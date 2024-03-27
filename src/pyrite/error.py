import sys

class PyriteError(Exception): 

    def __init__(self, message: str) -> None:
        self.message = message

class InputError(PyriteError): ...


def excepthook(exception_type, exception, traceback):

    if isinstance(exception, PyriteError):
        print(f"\nPyrite {exception_type.__name__}: {exception.message}")

    else:
        sys.__excepthook__(exception_type, exception, traceback)


sys.excepthook = excepthook