//! Rust -> Python logging bridge.
//!
//! This module forwards `tracing` events emitted anywhere in the native
//! core into the three Python loggers defined in `crisprme2/logger.py`:
//!
//! | tracing level   | Python destination            | file(s)                |
//! |-----------------|-------------------------------|------------------------|
//! | `ERROR`         | `errorlog`   (`.error`)       | `errors.log`           |
//! | `WARN`          | `basiclog` + `verboselog`     | `basic.log` + verbose  |
//! |                 | (`.info`, `"WARNING: "` prefix)|                       |
//! | `INFO`          | `basiclog` + `verboselog`     | `basic.log` + verbose  |
//! | `DEBUG`/`TRACE` | `verboselog` (`.debug`)       | `verbose.log`          |
//!
//! Design notes
//! ------------
//! * The performance-critical / domain code (e.g. `crispr::pam`) never
//!   touches Python, it only emits `tracing::{debug,info,error}!`. This
//!   layer is the single point of coupling to the interpreter.
//! * `WARN` is remapped onto `INFO` (with a textual marker) because the
//!   Python handlers use *exact-level* filters and expose no WARNING sink.
//! * Forwarding acquires the GIL, so it is intended for cold / low-rate
//!   paths (PAM parsing, pipeline milestones). The installed level filter
//!   drops `TRACE`, keeping hot-loop `trace!` events out of the interpreter.

use std::fmt::{Debug, Write as _};

use pyo3::prelude::*;
use pyo3::types::PyAnyMethods;
use tracing::field::{Field, Visit};
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::layer::{Context, Layer};
use tracing_subscriber::registry::LookupSpan;

/// Owned handles to the three underlying `logging.Logger` objects.
///
/// `Py<PyAny>` is `Send + Sync`, so the layer is safe to install as a
/// global subscriber shared across worker threads.
struct PyLoggers {
    basic: Py<PyAny>,   // -> basic.log   (INFO only)
    verbose: Py<PyAny>, // -> verbose.log (DEBUG + INFO)
    error: Py<PyAny>,   // -> errors.log  (ERROR)
}

/// A `tracing` [`Layer`] that forwards events to the Python loggers.
pub struct PyLoggerLayer {
    loggers: PyLoggers,
}

impl PyLoggerLayer {
    /// Build a layer from a Python `CrisprmeLoggers` bundle.
    ///
    /// Extracts the *underlying* `logging.Logger` from each wrapper via
    /// `get_logger()`, so forwarding writes records directly and respects
    /// the per-handler level filters configured in Python.
    ///
    /// # Errors
    /// Propagates any attribute/method error from the Python side (e.g. a
    /// bundle that is not a `CrisprmeLoggers`).
    pub fn from_bundle(bundle: &Bound<'_, PyAny>) -> PyResult<Self> {
        let basic = bundle.getattr("basiclog")?.call_method0("get_logger")?.unbind();
        let verbose = bundle.getattr("verboselog")?.call_method0("get_logger")?.unbind();
        let error = bundle.getattr("errorlog")?.call_method0("get_logger")?.unbind();
        Ok(Self { loggers: PyLoggers { basic, verbose, error } })
    }
}

/// Visitor that pulls the formatted `message` field out of an event.
#[derive(Default)]
struct MessageVisitor {
    message: String,
}

impl Visit for MessageVisitor {
    fn record_debug(&mut self, field: &Field, value: &dyn Debug) {
        if field.name() == "message" {
            // `tracing`'s message is recorded through the Debug channel.
            let _ = write!(self.message, "{value:?}");
        }
    }
}

impl<S> Layer<S> for PyLoggerLayer
where
    S: Subscriber + for<'a> LookupSpan<'a>,
{
    fn on_event(&self, event: &Event<'_>, _ctx: Context<'_, S>) {
        let mut visitor = MessageVisitor::default();
        event.record(&mut visitor);
        let msg = visitor.message;
        let level = *event.metadata().level();

        // Acquire the GIL and dispatch. Errors from the Python call are
        // intentionally swallowed: a logging failure must never unwind
        // through instrumentation.
        Python::attach(|py| {
            let emit = |logger: &Py<PyAny>, method: &str, text: &str| {
                let _ = logger.bind(py).call_method1(method, (text,));
            };
            match level {
                Level::ERROR => emit(&self.loggers.error, "error", &msg),
                Level::WARN => {
                    let warn = format!("WARNING: {msg}");
                    emit(&self.loggers.basic, "info", &warn);
                    emit(&self.loggers.verbose, "info", &warn);
                }
                Level::INFO => {
                    emit(&self.loggers.basic, "info", &msg);
                }
                Level::DEBUG | Level::TRACE => {
                    emit(&self.loggers.verbose, "debug", &msg);
                }
            }
        });
    }
}
