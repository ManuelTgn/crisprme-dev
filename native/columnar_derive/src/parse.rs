use syn::{Data, DeriveInput, Fields, Ident, LitStr, Type, spanned::Spanned};

/// Raw parsed information for a single struct field.
pub struct ParsedField<'a> {
    /// The field's identifier (e.g. `id`, `score`).
    pub name: &'a Ident,
    /// The field's declared type, as-is from the source.
    pub ty: &'a Type,
    /// `true` when `#[columnar(group)]` is present — the field will be
    /// expanded into N separate sub-columns.
    pub group: bool,
}

/// Raw parsed information for the whole struct being derived.
pub struct ParsedStruct<'a> {
    /// The struct's identifier (e.g. `Foo`).
    pub struct_name: &'a Ident,
    /// All named fields, in declaration order.
    pub fields: Vec<ParsedField<'a>>,
}

/// Extract field and struct attributes from a [`DeriveInput`].
/// Returns an error if the input is not a struct with named fields.
pub fn parse(input: &DeriveInput) -> syn::Result<ParsedStruct<'_>> {
    let fields = match &input.data {
        Data::Struct(s) => match &s.fields {
            Fields::Named(f) => &f.named,
            _ => {
                return Err(syn::Error::new(
                    input.span(),
                    "Columnar requires a struct with named fields",
                ));
            }
        },
        _ => {
            return Err(syn::Error::new(
                input.span(),
                "Columnar can only be derived for structs",
            ));
        }
    };

    let mut parsed_fields = Vec::new();
    for field in fields.iter() {
        let name = field.ident.as_ref().unwrap();
        let ty = &field.ty;

        let mut group = false;
        for attr in &field.attrs {
            if attr.path().is_ident("columnar") {
                let _ = attr.parse_nested_meta(|meta| {
                    if meta.path.is_ident("group") {
                        group = true;
                    }
                    Ok(())
                });
            }
        }

        parsed_fields.push(ParsedField { name, ty, group });
    }

    Ok(ParsedStruct {
        struct_name: &input.ident,
        fields: parsed_fields,
    })
}
