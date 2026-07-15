use proc_macro2::Span;
use syn::{Ident, Type};

use crate::parse::ParsedStruct;

/// Intermediate representation for a single field.
/// Both scalar, array and `#[columnar(group)]` fields share this type.
pub struct FieldIR<'a> {
    /// The field's identifier.
    pub name: &'a Ident,
    /// The per-element type: `T` from `[T; N]` for arrays and groups,
    /// or the field type itself for scalar and array fields.
    pub ty: &'a Type,
    /// How to manage this field: scalar, array or group
    pub kind: FieldKindIR<'a>,
}

pub enum FieldKindIR<'a> {
    Scalar,
    Array { len: &'a syn::Expr },
    Group { len: &'a syn::Expr },
}

/// Intermediate representation for the whole schema.
pub struct SchemaIR<'a> {
    /// The original struct's identifier (e.g. `Foo`).
    pub struct_name: &'a Ident,
    /// All fields, in declaration order.
    pub fields: Vec<FieldIR<'a>>,
}

/// Lower a [`ParsedStruct`] into a [`SchemaIR`].

pub fn lower<'a>(parsed: &ParsedStruct<'a>) -> syn::Result<SchemaIR<'a>> {
    let mut fields = Vec::new();
    for f in &parsed.fields {
        match f.ty {
            Type::Array(arr) => {
                let kind = match f.group {
                    true => FieldKindIR::Group { len: &arr.len },
                    false => FieldKindIR::Array { len: &arr.len },
                };
                fields.push(FieldIR {
                    name: f.name,
                    ty: &*arr.elem,
                    kind,
                });
            }
            _ => {
                fields.push(FieldIR {
                    kind: FieldKindIR::Scalar,
                    name: f.name,
                    ty: f.ty,
                });
            }
        }
    }

    Ok(SchemaIR {
        struct_name: parsed.struct_name,
        fields,
    })
}
