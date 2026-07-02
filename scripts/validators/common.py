from functools import wraps
from pathlib import Path
from typing import Any, Callable, Protocol, cast

from tqdm import tqdm

from finds.constants import VALIDATION_OUTPUT_PATH


class ValidationFunction(Protocol):
    """
    A custom protocol for any validation function. Allows for the addition of
    :code:`name` and :code:`filename` attributes.

    Attributes
    ==========
    name: str
      The name of the validation, like a desired figure name for reproduction.

    filename: Path
      The output filename for the produced figure, video, etc.
    """
    name: str
    filename: Path

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        pass


def produces_validation(*, name: str, output_type: str = 'png'):
    """
    Simple internal decorator for all validation functions that produce
    figures.

    :param name: The name of the validation.
    :type  name: str

    :param output_type: The output type for the produced figure, video, etc.
      (default is :code:`.png`).
    :type  output_type: str
    """
    def decorator(function_raw: Callable[..., Any]):
        function = cast(ValidationFunction, function_raw)

        function.name = name
        function.filename = VALIDATION_OUTPUT_PATH /\
            f'validation-{name}.{output_type}'

        @wraps(function)
        def wrapper(*args, **kwargs):
            tqdm.write(f"Beginning validation of figure {function.name}")
            result = function(*args, **kwargs)
            tqdm.write(f"Generated validation to {function.filename}")
            return result

        return wrapper

    return decorator
