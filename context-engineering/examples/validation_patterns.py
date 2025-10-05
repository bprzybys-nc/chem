"""
Validation Pattern Examples

This file demonstrates common validation patterns used in the project.
"""

def validate_with_fast_failure(data: dict) -> dict:
    """
    Fast-failure validation pattern.

    Validates data and raises specific exceptions immediately on error.
    No silent failures or default values.
    """
    if not data:
        raise ValueError("Data cannot be empty")

    if "required_field" not in data:
        raise KeyError("Missing required_field. Available fields: " + ", ".join(data.keys()))

    return data


def validate_with_type_checking(value: str, expected_type: type) -> bool:
    """
    Type checking validation pattern.

    Validates value type and provides actionable error messages.
    """
    if not isinstance(value, expected_type):
        raise TypeError(
            f"Expected {expected_type.__name__}, got {type(value).__name__}. "
            f"Please convert value to {expected_type.__name__} before passing."
        )

    return True
