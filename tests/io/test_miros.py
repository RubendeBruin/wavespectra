import os
import xarray as xr
import shutil
import pytest
from tempfile import mkdtemp

from wavespectra import read_miros
from wavespectra.core.attributes import attrs

FILES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sample_files")


class TestMiros(object):
    """Test Octopus writer."""

    @classmethod
    def setup_class(self):
        """Setup class."""
        self.filename = os.path.join(FILES_DIR, "MIR_WaveDta_20210502_DAY.DF038")


    def test_read_miros(self):
        dset = read_miros(self.filename)
