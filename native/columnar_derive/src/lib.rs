use proc_macro::TokenStream;
use syn::{DeriveInput, parse_macro_input};

mod codegen;
mod parse;
mod repr;

#[proc_macro_derive(Columnar, attributes(columnar))]
pub fn derive_columnar_schema(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    match expand(&input) {
        Ok(ts) => ts.into(),
        Err(e) => e.to_compile_error().into(),
    }
}

fn expand(input: &DeriveInput) -> syn::Result<proc_macro2::TokenStream> {
    let parsed = parse::parse(input)?;
    let ir = repr::lower(&parsed)?;
    Ok(codegen::generate(&ir))
}
