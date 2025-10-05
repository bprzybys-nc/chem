"""
Error Handling Pattern Examples

This file demonstrates proper error handling patterns.
"""

def process_with_explicit_errors(data: dict) -> dict:
    """
    Explicit error handling pattern.

    - No fishy fallbacks
    - No silent failures
    - Actionable error messages
    - Fast failure on invalid data
    """
    try:
        result = perform_processing(data)
        return result
    except KeyError as e:
        raise ValueError(
            f"Processing failed - missing key: {e}. "
            f"Required keys: field1, field2, field3"
        )
    except TypeError as e:
        raise TypeError(
            f"Processing failed - invalid type: {e}. "
            f"Expected dict with string values."
        )


def process_with_troubleshooting(data: dict) -> dict:
    """
    Error handling with troubleshooting guidance.

    Provides clear troubleshooting steps in error messages.
    """
    try:
        result = perform_processing(data)
        return result
    except Exception as e:
        error_msg = (
            f"Processing failed: {e}\n"
            f"Troubleshooting steps:\n"
            f"  1. Check if data format matches expected schema\n"
            f"  2. Verify all required fields are present\n"
            f"  3. Ensure field values are correct types"
        )
        raise RuntimeError(error_msg)


# FIXME: Placeholder function for demonstration
def perform_processing(data: dict) -> dict:
    """Placeholder processing function."""
    return data
