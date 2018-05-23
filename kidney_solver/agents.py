"""Various agents in the KE model."""


class UnknownAgeException(Exception):
    """Raised if the age of an agent is not known."""
    pass


class Donor(object):
    """A donor in a KE."""

    def __init__(self, ident, altruistic, age=None):
        self._ident = ident
        self._altruistic = altruistic
        self._age = age

    def age(self):
        """The age of the donor.

        Throws an exception if the age is not known.
        """
        if self._age is None:
            raise UnknownAgeException()
        return self._age

    def index(self):
        """The index of the donor.
        """
        return self._ident

    def is_altruistic(self):
        """Is this donor altruistic?"""
        return self._altruistic

    def __repr__(self):
        return str(self)

    def __str__(self):
        if self._altruistic:
            return "D(A)%s" % self.index()
        else:
            return "D%s" % self.index()
