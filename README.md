## LCMAP Change Detection

This project contains application source code for Change Detection C library
and related scripts.


## Installation

To install, simply run the top-level ``make`` target:

```bash
$ git clone git@github.com:USGS-EROS/lcmap-change-detection-c.git
$ cd lcmap-change-detection-c
$ make

ccdc depends on the GNU Scientific Library (GSL)
(https://www.gnu.org/software/gsl/).
These libs are used to replace some of the mathematical funcions in the
proprietary Matlab (mathworks.com) version.  A local copy of multirobust.c
has been modified and used so other applications can link to a clean copy of
the gsl lib.
The Matalb-formatted output file(s) (.mat) depend on the MAT File I/O Library:
https://sourceforge.net/projects/matio/files/matio/1.5.2/

Therefore, in the ccdc source sub-directory, make.sh has been added to define
these (and any other dependencies) using environment variables in the call to
make, which uses Makefile.
```

The executables and scripts will be installed into ``./bin`` by default. This
can be overridden by setting a ``BIN`` environment variable or using a ``BIN``
variable when running the target:

```bash
$ BIN=/my/path/bin make
```

For additional notes, such as installing dependencies (Ubuntu), overriding
``Makefile`` variables, etc., see:

* [Building CCDC](../..//wiki/Building-CCDC)


## Usage

[We're in active development on making this not only work, but be usable.
Ticket #5 has some early usability notes + tasks that we're trying to hit right
away, if you're interested in tracking this.]

Continuous Change Detection and Classification
Version 05.04, topic/tile-based

usage:<br>
ccdc --row=<input row number> --col=<input col number> --num-rows=<number of rows> --num-cols=<number of columns> [--in-path=<input directory>] [--out-path=<output directory>] [--data-type=tifs|bip] [--scene-list-file=<file with list of sceneIDs>] [--verbose]

where the following parameters are required:<br>
    --row=: input row number<br>
    --col=: input col number

and the following parameters are optional:<br>
    --num-rows=: input number of rows<br>
    --num-cols=: input number of columns<br>
    --in-path=: input data directory location<br>
    --out-path=: directory location for output files<br>
    --data-type=: type of input data files to ingest<br>
    --scene-list-file=: file name containing list of sceneIDs (default is all files in in-path)<br>
    -verbose: should intermediate messages be printed? (default is false)

ccdc --help will print the usage message


Example:<br>
ccdc --row=3845 --col=2918 --num-rows=1 --num-cols=1 --in-path=/data/user/in --out-path=/home/user/out --data-type=bip --scene-list-file=/home/user/scene_list.txt --verbose

An example of how to pipe input from stdin and output to stdout:<br>
ccdc --row=3845 --col=2918 --in-path=stdin --out-path=stdout --verbose < pixel_value_text_file.txt > coeffs_results_text_file.txt

The stdout option eliminates the creation of the output binary file,<br>
coeffs are just printed to stdout.  This can be re-directed to a text file,<br>
or piped to another program.

Note: Previously, the ccdc had to be run from the directory where the input data<br>
      are located.  Now, input and output directory location specifications are<br>
      used.  If in-path or out-path are not specified, current working directory<br>
      is assumed.  If scene-file-name is not specified, all scenes in in-path<br>
      are processed.


## Input Data

The filesystem data located in in-path is expected to be gridded Landsat<br>
Analysis Ready Data (ARD), co-registered and radiometrically correlated.<br>
For landsat 8, Surface reflectance (SR) bands 2-7 are used for the spectral<br>
image bands, and Top Of Atmosphere (TOA) band 10 is use for the thermal band,<br>
in addition to the Cloud Mask band (cfmask). for Landsats 4, 5, and 7 SR bands<br>
1-5 and 7 are used for spectral imagery, and TOA band 6 for thermal data.<br>

For --data-type=bip, these bands are stacked into a single<br>
Band Interleaved by Pixel (BIP) format file, with an ENVI style header file<br>
with the same name.  The base of the file name is the sceneID.  There is<br>
a pair of these files for each sceneID defined in the input specifications.<br>
This data format is compatible with the original matalb-developed application.<br>

For example:

LC80440272013106/LC80440272013106LGN01_MTLstack<br>
LC80440272013106/LC80440272013106LGN01_MTLstack.hdr<br>
.......


For --data-type=tifs, each band is in a separate .tif file, one set of band<br>
files for each sceneID:

LC80440272013106LGN01_cfmask.img<br>
LC80440272013106LGN01_sr_band2.hdr<br>
LC80440272013106LGN01_sr_band2.img<br>
LC80440272013106LGN01_sr_band3.img<br>
LC80440272013106LGN01_sr_band4.img<br>
LC80440272013106LGN01_sr_band5.img<br>
LC80440272013106LGN01_sr_band6.img<br>
LC80440272013106LGN01_sr_band7.img<br>
LC80440272013106LGN01_toa_band10.img<br>
.......


The values required for --in-path=stdin are a set values per line for each<br>
input sceneID: Julian date (since 0); the 6 spectral band values; thermal band<br>
value; and cfmask value.  For example:<br>
2456397 3737 3578 3524 3833 1343 1377 2502 4<br>
2456413 596 808 898 1408 220 173 2696 4<br>
2456429 7745 7945 8147 8537 5463 4454 2649 4<br>
2456445 94 156 164 758 807 492 2809 0<br>
2456461 264 564 467 1680 1465 901 2866 0<br>
2456477 150 299 245 1308 1108 650 2950 0<br>
2456493 156 297 231 1377 1085 629 2965 0<br>
2456509 2792 2171 2242 3050 1385 1158 2683 4<br>
2456525 115 220 180 1014 809 465 2886 0<br>
2456541 4732 4605 4570 4892 2174 2287 2349 4<br>
2456557 73 85 65 410 208 100 2800 2<br>
2456573 71 121 106 610 417 278 2774 0<br>
2456589 54 81 64 229 170 101 2800 1<br>
......

## Development

Development notes for C-CCDC are maintained in the project wiki. For more
details, see:

 * [CCDC Development](../../wiki/CCDC Development)
 * [Using CCDC with Docker](../../CCDC-%26-Docker)


## Implementation

### CCDC - Continuous Change Detection and Classification (Algorithm)

* <b>NOTE:</b> This algorithm is not validated and considered prototype.
* See [CCDC ADD](http://landsat.usgs.gov/documents/ccdc_add.pdf) for the
  detailed description.


## More Information

This project is hosted by the US Geological Survey (USGS) Earth Resources
Observation and Science (EROS) Land Change Monitoring, Assessment, and
Projection ([LCMAP](https://github.com/USGS-EROS?utf8=%E2%9C%93&query=lcmap))
Project.  For questions regarding this source code, please contact the
[Landsat Contact Us](https://landsat.usgs.gov/contactus.php) page and specify
``USGS LCMAP`` in the "Regarding" section.
