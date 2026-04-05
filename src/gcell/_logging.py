import contextlib
import logging
import sys

from . import GRADIO_AVAILABLE

with contextlib.suppress(ImportError):
    import gradio as gr


class GradioLoggerAdapter(logging.LoggerAdapter):
    """
    A logger adapter that checks if code is running in a Gradio context
    and raises Gradio errors when appropriate.
    """

    def __init__(self, logger: logging.Logger, extra: dict | None = None):
        super().__init__(logger, extra or {})

    def error(self, msg: str, *args, **kwargs):
        """
        Log an error message. If running in a Gradio context, raise a Gradio error.
        Otherwise, log the error and raise a ValueError.

        Args:
            msg: The error message to log
            *args: Additional positional arguments for the logger
            **kwargs: Additional keyword arguments for the logger

        Raises:
            gr.Error: If running in a Gradio context
            ValueError: If not running in a Gradio context
        """
        self.logger.error(msg, *args, **kwargs)
        if GRADIO_AVAILABLE and hasattr(gr, 'get_state') and gr.get_state() is not None:
            raise gr.Error(msg)
        raise ValueError(msg)


def get_logger(name: str) -> GradioLoggerAdapter:
    """
    Get a logger instance with the given name that's aware of Gradio context.

    Args:
        name: The name of the logger

    Returns:
        GradioLoggerAdapter: A logger adapter that handles both Gradio and non-Gradio contexts
    """
    logger = logging.getLogger(name)

    # Only add handlers if none exist
    if not logger.handlers:
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(console_handler)

        # Set default level to INFO
        logger.setLevel(logging.INFO)

    return GradioLoggerAdapter(logger)
