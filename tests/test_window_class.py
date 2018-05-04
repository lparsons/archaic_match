import pytest
from archaic_match.classmodule import Window


@pytest.fixture
def w():
    return Window("chr", 0, 10)


def test_create_window(w):
    assert isinstance(w, Window)


def test_window_size(w):
    assert w.size == 10


def test_window_size_0(w):
    with pytest.raises(AttributeError):
        Window("chr", 10, 10)


def test_window_negative_start(w):
    with pytest.raises(AttributeError):
        Window("chr", -1, 10)


def test_window_isc_is_none(w):
    assert w.informative_sites_count is None


def test_window_isf_is_none(w):
    assert w.informative_sites_frequency is None


def test_window_set_isc(w):
    w.informative_sites_count = 1
    assert w.informative_sites_count == 1


def test_window_set_isc_too_large(w):
    with pytest.raises(AttributeError):
        w.informative_sites_count = 11


def test_window_set_isf(w):
    w.informative_sites_count = 1
    assert w.informative_sites_frequency == 0.1
