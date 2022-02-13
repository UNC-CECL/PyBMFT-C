import pkg_resources

from .bmftc import Bmftc

__version__ = pkg_resources.get_distribution("bmftc").version
__all__ = ["Bmftc"]

del pkg_resources