## LCMAP Change Detection

This project contains application source code for Change Detection C library
and related scripts.


## Installation

To install, simply run the top-level ``make`` target:

```bash
$ git clone git@github.com:USGS-EROS/lcmap-change-detection-c.git
$ cd lcmap-change-detection-c
$ make
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
