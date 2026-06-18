from functools import wraps
from pathlib import Path
from typing import Callable

from tqdm import tqdm

#: The validation output directory
VALIDATION_OUTPUT_DIR = Path(f'output/validation')


def produces_validation(*, name: str, output_type: str = 'png'):
    """
    Simple internal decorator for all validation functions that produce
    figures.
    """
    def decorator(function: Callable) -> Callable:
        function.name = name
        function.filename = VALIDATION_OUTPUT_DIR /\
            f'validation-{name}.{output_type}'

        @wraps(function)
        def wrapper(*args, **kwargs):
            tqdm.write(f"Beginning validation of figure {function.name}")
            result = function(*args, **kwargs)
            tqdm.write(f"Generated validation to {function.filename}")
            return result

        return wrapper

    return decorator
