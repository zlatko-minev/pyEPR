"""
Tests for DistributedAnalysis._variation_index — the variation-label-to-integer
coercion that replaces the old ureg(variation) approach.

ureg(variation) was broken because:
  - On some pint versions it returns a Quantity instead of an int, causing
    TypeError when used to index a tuple.
  - It cannot handle full descriptor strings like "Cj='2fF' Lj='12nH'"
    (raises pint.UndefinedUnitError).

These tests verify the replacement logic in isolation (no HFSS session needed).
"""
import types
import pytest

from pyEPR.core_distributed_analysis import DistributedAnalysis


def _make_stub(list_variations, variations=None):
    """Return a minimal object that satisfies _variation_index's attribute needs."""
    obj = types.SimpleNamespace()
    obj._list_variations = list_variations
    obj.variations = variations or [str(i) for i in range(len(list_variations))]
    obj._variation_index = DistributedAnalysis._variation_index.__get__(obj)
    return obj


class TestVariationIndex:

    def setup_method(self):
        descriptors = (
            "Cj='2fF' Lj='12nH'",
            "Cj='2fF' Lj='12.5nH'",
            "Cj='2fF' Lj='13nH'",
        )
        self.stub = _make_stub(descriptors)

    # --- fast path: digit strings and ints ---

    def test_digit_string(self):
        assert self.stub._variation_index('0') == 0
        assert self.stub._variation_index('1') == 1
        assert self.stub._variation_index('2') == 2

    def test_integer(self):
        assert self.stub._variation_index(0) == 0
        assert self.stub._variation_index(2) == 2

    # --- fallback: full descriptor strings ---

    def test_full_descriptor_string(self):
        assert self.stub._variation_index("Cj='2fF' Lj='12nH'") == 0
        assert self.stub._variation_index("Cj='2fF' Lj='13nH'") == 2

    # --- fallback: variation label from self.variations ---

    def test_variation_label_fallback(self):
        # self.variations = ['0', '1', '2'] — same as digit strings in this case,
        # but the fallback path is exercised when _list_variations is a non-tuple
        # that doesn't contain the label.
        stub = _make_stub(("a", "b", "c"), variations=["x0", "x1", "x2"])
        assert stub._variation_index("x0") == 0
        assert stub._variation_index("x2") == 2

    # --- error path ---

    def test_unknown_variation_raises(self):
        with pytest.raises(ValueError, match="variation="):
            self.stub._variation_index("not-a-valid-variation")

    def test_none_raises(self):
        # _variation_index should never be called with None (callers guard it),
        # but if it is, it should raise rather than silently return 0.
        with pytest.raises((TypeError, ValueError)):
            self.stub._variation_index(None)

    # --- regression: the ureg bug ---

    def test_pint_ureg_would_have_failed_on_descriptor(self):
        """Demonstrate that ureg cannot parse a descriptor string."""
        from pyEPR.ansys import ureg
        with pytest.raises(Exception):
            _ = ureg("Cj='2fF' Lj='12nH'")

    def test_pint_ureg_digit_string_instability(self):
        """ureg('0') may return Quantity(0) not int — verify _variation_index is stable."""
        from pyEPR.ansys import ureg
        result = self.stub._variation_index('0')
        assert result == 0
        assert type(result) is int
