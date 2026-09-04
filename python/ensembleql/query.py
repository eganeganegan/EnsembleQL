"""Query validation helpers."""


def validate(text: str) -> None:
    """Raise QueryError if *text* is not valid EnsembleQL syntax."""
    try:
        from ._core import parse_query
    except ImportError as exc:
        raise RuntimeError("EnsembleQL's native extension is not built") from exc
    parse_query(text)
