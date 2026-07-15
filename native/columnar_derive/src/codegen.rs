use proc_macro2::{Span, TokenStream};
use quote::{format_ident, quote};
use syn::Ident;

use crate::repr::{FieldKindIR, SchemaIR};

/// Top-level code generation entry point.
pub fn generate(ir: &SchemaIR) -> TokenStream {
    // Number of slots required for this schema
    let total_slots = compute_total_slots(ir);

    let schema = generate_schema(ir, total_slots);
    let cols = generate_cols(ir, total_slots);

    quote! {
        #schema
        #cols
    }
}

/// Count how many slots a field uses.
fn field_slot_count(kind: &FieldKindIR) -> usize {
    match kind {
        FieldKindIR::Scalar | FieldKindIR::Array { .. } => 1,
        FieldKindIR::Group { len } => {
            let tokens = quote! { #len };
            syn::parse2::<syn::LitInt>(tokens)
                .expect("Group length must be an integer literal")
                .base10_parse::<usize>()
                .expect("Group length must be a valid usize")
        }
    }
}

fn compute_total_slots(ir: &SchemaIR) -> usize {
    ir.fields.iter().map(|f| field_slot_count(&f.kind)).sum()
}

fn generate_schema(ir: &SchemaIR, total_slots: usize) -> TokenStream {
    let struct_name = ir.struct_name;
    quote! {
        unsafe impl ::columnar::Schema for #struct_name {
            const SLOTS: usize = #total_slots;
        }
    }
}

fn generate_cols(ir: &SchemaIR, total_slots: usize) -> TokenStream {
    let struct_name = ir.struct_name;
    let cols_name = Ident::new(&format!("{}Cols", ir.struct_name), Span::call_site());

    let typed_struct_name = Ident::new(&format!("{}Frame", ir.struct_name), Span::call_site());

    // Compute cumulative slot offsets
    let mut slot_offsets: Vec<usize> = Vec::new();
    let mut offset = 0usize;
    for f in &ir.fields {
        slot_offsets.push(offset);
        offset += field_slot_count(&f.kind);
    }

    // Generate the Cols struct fields
    let cols_struct_fields: Vec<_> = ir
        .fields
        .iter()
        .map(|f| {
            let name = f.name;
            let ty = f.ty;
            match &f.kind {
                FieldKindIR::Scalar => {
                    quote! { pub #name: ::columnar::Column<'a, #ty> }
                }
                FieldKindIR::Array { len } => {
                    quote! { pub #name: ::columnar::Column<'a, [#ty; #len]> }
                }
                FieldKindIR::Group { len } => {
                    quote! { pub #name: ::columnar::ColumnGroup<'a, #ty, #len> }
                }
            }
        })
        .collect();

    // Generate field initializers using unsafe raw pointer indexing
    let field_inits: Vec<_> = ir
        .fields
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let name = f.name;
            let off = slot_offsets[i];
            match &f.kind {
                FieldKindIR::Scalar | FieldKindIR::Array { .. } => {
                    quote! { #name: ::columnar::Column::new(unsafe { &mut *__ptr.add(#off) }) }
                }
                FieldKindIR::Group { .. } => {
                    let n = field_slot_count(&f.kind);
                    let group_slots: Vec<_> = (0..n)
                        .map(|j| {
                            let idx = off + j;
                            quote! { unsafe { &mut *__ptr.add(#idx) } }
                        })
                        .collect();
                    quote! { #name: ::columnar::ColumnGroup::new([#( #group_slots ),*]) }
                }
            }
        })
        .collect();

    // Generate alloc calls for `new()`
    let alloc_calls: Vec<_> = ir
        .fields
        .iter()
        .map(|f| {
            let name = f.name;
            quote! { __cols.#name.alloc(__pool, __rows); }
        })
        .collect();

    quote! {

        pub struct #cols_name<'a> {
            #( #cols_struct_fields, )*
        }

        /// Wrap the TypedFrame into a custom struct for foreign implementation
        pub struct #typed_struct_name(::columnar::TypedFrame<#struct_name>);
        impl #typed_struct_name {

            /// Create an empty frame with this schema
            pub fn empty() -> Self {
                Self(::columnar::TypedFrame {
                    _schema: ::std::marker::PhantomData,
                    frame: ::columnar::frame::DynFrame::empty(
                        #struct_name::SLOTS),
                })
            }

            /// Allocate frame for `rows` rows
            pub fn alloc(__pool: &::columnar::MemoryPool, __rows: usize) -> Self {
                let mut __result = Self::empty();
                __result.with_cols(|mut __cols| {
                    #( #alloc_calls )*
                });
                __result
            }

            /// Mutable access to the inner frame (used by derive macro)
            pub fn frame_mut(&mut self) -> &mut ::columnar::frame::DynFrame { &mut self.0.frame }

            /// Attach a schema to an existing frame.
            /// # Safety
            /// The caller must ensure the frame's layout matches `S`.
            pub unsafe fn attach(frame: ::columnar::frame::DynFrame) -> Self {
                Self(::columnar::TypedFrame {
                    _schema: ::std::marker::PhantomData,
                    frame
                })
            }

            /// Consume the wrapper and return the inner frame
            pub fn detach(self) -> ::columnar::frame::DynFrame { self.0.frame }

            /// Expose columns of typed frame
            pub fn with_cols<F, R>(&mut self, f: F) -> R
            where
                F: FnOnce(#cols_name<'_>) -> R
            {
                let __ptr = unsafe { self.frame_mut().slots_ptr() };
                let __cols = #cols_name {
                    #( #field_inits, )*
                };
                f(__cols)
            }
        }

        impl ::columnar::Share for #typed_struct_name {
            fn share(&mut self) -> Self {
                Self(self.0.share())
            }
        }
    }
}
