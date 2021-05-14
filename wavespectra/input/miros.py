import xarray as xr
import pandas as pd
import numpy as np

from wavespectra.specdataset import SpecDataset

def read_miros(filename):
    """
    Miros standard file name convention for DF-038 history files is:
    File
    length
    File name convention
    HOUR SSS_XXX_YYYYMMDD_hhmm_HRS.DF038
    DAY SSS_XXX_YYYYMMDD_DAY.DF038
    MONTH SSS_XXX_YYYYMM_MON.DF038
    Note! The minute field (mm) of the time stamp in an hour file name should be 00, i.e. the
    data history should only include files with names referring to full hours.
    Here:
    SSS = Site Identifier (three alpha-numeric characters, letters are written in upper case).
    XXX = Sequence Number (free text limited to letters and digits).
    YYYY = Year (four digits).
    MM = Month (two digits with a leading zero if required).
    DD = Day (two digits with a leading zero if required).
    hh = Hours (two digits in 24 hours format with a leading zero if required).
    mm = Minutes (two digits with a leading zero if required).
    HRS/
    DAY/
    MON
    = File length:
     HRS for files containing data from a period of up to one hour.
     DAY for files containing data from a period of up to one day (24 hours).
     MON for files containing data from a period of up to one month (up to
    31 days.
    DF038 = File Extension (letters are written in upper case).


    Directional information is found in the Parameter identification section, through the four
    parameters:
     Direction_Convention
     Direction_Relative_To
     Direction_Offset
     Number_Of_Directions
    and in the Heading column of the Data section.


    Direction_Convention=Approaching
    Direction_Relative_To=Heading
    Direction_Offset=277
    Number_Of_Directions=036
    Number_Of_Frequencies=037
    Start_Frequency=0.0312500
    Frequency_Resolution=0.0078125
    Number_Of_Status_Codes=1
    Heading=Yes
    Position=Yes

    """

    header_rows = 0
    metadata = {}
    with open(filename, 'r') as file:
        # recall that with ... as ... opens (and closes!) the file

        while header_rows < 30:
            line = file.readline()
            line = line.strip('\n')  # strip the eol
            # print(line)

            header_rows += 1

            # do we have the start of the data?
            if '[DATA]' in line:
                break

            # can we break the line into a x=y ?
            parts = line.split('=')
            if len(parts) == 2:
                metadata[parts[0]] = parts[1]

    n_directions = int(metadata['Number_Of_Directions'])
    n_frequencies = int(metadata['Number_Of_Frequencies'])

    start_freq = float(metadata['Start_Frequency'])
    freq_resolution = float(metadata['Frequency_Resolution'])

    freqs = [(start_freq + i * freq_resolution) for i in range(0, n_frequencies)]

    # Read the [DATA] section
    #
    # Parse the first three columns as date (so including the 'Z' for time-zone )
    data = pd.read_csv(filename, header=None, sep='\s+', skiprows=header_rows,
                       index_col=0, parse_dates=[[0, 1, 2]])

    # get the time
    times = data.index.values

    # get the values as np array
    values = data.values

    # trim the data to the directional spectrum
    #
    # trim the point-spectrum from the start
    i_start = n_frequencies
    i_end = i_start + n_directions * n_frequencies

    dir_spectra = values[:, i_start:i_end]

    # loop over the rows
    spectra = []
    for spectrum in dir_spectra:
        spectra.append(spectrum.reshape(n_directions, n_frequencies))

    # Now read the remaining columns one-by-one
    #
    # Depending on the metadata, some data may not be present

    icol = i_end

    # Status column
    # Note that
    # Number_Of_Status_Codes can be set to 000, implying that the status
    # column is empty.
    if int(metadata['Number_Of_Status_Codes']) > 0:
        icol += 1  # we do not read the status-codes

    if 'No' in metadata['Heading']:
        heading = None  # no heading given
    else:
        heading = values[:, icol]
        icol += 1
        heading_spread = values[:, icol]
        icol += 1

    if 'No' in metadata['Position']:
        latt = None
        long = None
    else:
        latt = values[:, icol]
        icol += 1
        long = values[:, icol]
        icol += 1

    """
    Now determine the heading-grid


    4.2 Interpreting the directional information
    4.2.1 The Direction_Convention parameter
    The Direction_Convention parameter determines if the wave reference directions are given as
    Approaching or Receding, or “coming from” and “going to”, respectively. Switching between
    these two conventions is done by rotating the spectrum 180˚.

    """

    if 'pproaching' in metadata['Direction_Convention'] or \
        'coming from' in metadata['Direction_Convention']:
        convention = 0
    else:
        convention = 180

    """"
    4.2.2

    The number of directions determines the directional resolution.
    If Number_Of_Directions is set to 36, the directional resolution is 360/36 = 10˚. In this case the 36 direction columns of the Data
    section can be distributed as shown in figure 4.1

    Figure 4.1 shows direction 1 to be the sector between 355 and 5 degrees

    """

    heading_resolution = 360 / n_directions
    dirs = np.arange(0, 360, step=heading_resolution) + convention # note: excludes 360 itself

    """
    4.2.3 The Direction_offset parameter - Align the spectrum with the direction
    reference
    The Direction_Offset parameter indicates how much the spectrum must be rotated clockwise to
    be aligned with the direction reference, Heading or North, specified in the Direction_Relative_To
    parameter
    """

    try:
        direction_offset = float(metadata['Direction_Offset'])
    except:
        direction_offset = 0

    """
    4.2.4 Fixed heading - Direction_Relative_To=Heading
    This section describes the case when the source of the spectrum is a sensor mounted on a
    fixed site with a defined heading. As described in section 4.2.3 Direction_offset indicates how
    much the spectrum must be rotated clockwise to be aligned with the heading. The offset
    between heading and north is given in the heading column of the data section. Rotate the
    spectrum clockwise with this heading value to align it with north
    """

    if metadata['Direction_Relative_To'] == 'Heading':

        # now this is annoying. The heading varies which means we need to adjust the spectral information for each of the
        # spectra to align with the grid defined in dirs

        from scipy.interpolate import interp1d

        rotated_specs = []
        for head, spec in zip(heading, spectra):
            total = np.sum(spec)

            # spec ( heading , frequency )
            source_dir = np.mod(dirs + head, 360)
            source_dir = np.hstack([source_dir[-1] - 360, source_dir, source_dir[0] + 360])
            spec = np.vstack([spec[-1, :], spec, spec[0, :]])

            interpolant = interp1d(source_dir, spec, axis=0)
            rotated_spec = interpolant(dirs)

            # maintain the sum
            rotated_spec = rotated_spec * total / np.sum(rotated_spec)
            rotated_specs.append(rotated_spec)

    elif metadata['Direction_Relative_To'] == 'North':
        dirs = dirs + direction_offset

    dir_bin_width = heading_resolution  # convert to rad

    # convert S to s ?
    spectra = np.array(spectra, dtype=float)
    density = spectra  / dir_bin_width  # Correct to only divide by heading-resolution?

    from wavespectra.core.attributes import attrs

    coords = ((attrs.TIMENAME, times),
              (attrs.DIRNAME, dirs),
              (attrs.FREQNAME, freqs))
    dimensions = tuple([c[0] for c in coords])

    dset = xr.DataArray(density, coords=coords, dims=dimensions, name=attrs.SPECNAME
        ).to_dataset(name='efth')

    if latt is not None:
        dset[attrs.LATNAME] = ('time', latt)
    if long is not None:
        dset[attrs.LONNAME] = ('time', long)
    if heading is not None:
        dset['Miros_heading'] = ('time', heading)

    return dset

