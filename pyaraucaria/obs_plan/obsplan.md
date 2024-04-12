# Observing Plan Syntax

## Objectives

After evaluation of pros and cons, we decided to introduce new file format for observation plans, which will be used by the new (2022+) OCA observatory software.
As, in some simplification, observations plan is a sequence of command, the format of Observations Plan file is basically a syntax of kind of programming language.
because we want the new format to be human-readable, -editable and -manageable, we decide to introduce new language instead of using e.g. JSON syntax. 
New format should be somehow similar to ols obsplan files for the **IRIS** and **V16** telescopes and in general to be readable and understandable by astronomers without 
need to RTFM (*read the fantastic manual*).

We want to keep the syntax well-defined and intuitively but uniquely transformable to python dictionary object (and therefore to JSON, YAML etc.).
Moreover, we want to provide any developers with python package for parsing anf formatting those files.

## Introduction

We want the new language to act as internal language of the software, so the simplest form is just single command to be executed rather than full-fledged observation plan.
E.g. we want following text to be proper "program":
```
  WAIT ut=23:02:23
```
or, a bit more complicated single command
```
  OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=5/Ic/60,5/V/70
```
Both of single command lines are understandable for OCA astronomers 
(wait for UT 23:02:23, make photos of object HD193901 and ra/dec coordinates 20:23:35.8 -21:22:14.0 five times for
50 seconds with filter *Ic* 
and five times for 70 seconds with filter *V*). Syntactical break down of the last command goes as follows:
* command: `OBJECT`
* positional arguments (later Args): `HD193901`, `20:23:35.8`, `-21:22:14.0` - some from the tail may be optional 
* keyword arguments (later Kwargs): argument named `sequence` of value `5/I/60,5/V/70`. all keyword arguments can be optional

### Sequences
Sequence is a number of commands each written in unique, single line e.g:
```
WAIT ut=16:00:00
ZERO seq=15/Ic/0
DARK ZZ01 seq=10/V/300,10/Ic/200
DOMEFLAT AL2 seq=7/V/20,7/u/20
DOMEFLAT AL3 seq=10/i/100
SKYFLAT HD23 alt=60.0 az=270.0 seq=10/r/20,10/V/30 
SKYFLAT HD24 seq=10/g/a,10/V/a
WAIT wait=600
FOCUS NG31 12:12:12 20:20:20
OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/Ic/60,5/V/70
OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20
```


## Syntax Details

### Commands

#### OBJECT
* __*Description*__: Main observing plan command. If `ra`, `dec` or `alt`, `az` is not given, but `seq` is given, program
will do `seq` in actual telescope position. This command is possible to dithering in option.
* __*Default*__:  Coordinates epoch is J2000, dome is follow telescope, 
mirror covers are not opening/closing, dome cover is not opening/closing, dithering is OFF, tracking is ON.
* __*Args*__: `object_name`, `ra`, `dec`
* __*Kwargs*__: `seq`, `az`, `alt`, `tracking`, `dither`, `dome_follow`, `mirror_cover`, `epoch`, `observer`, `uobi`
* __*Syntax*__: `OBJECT [object_name] [optional: ra] [optional: dec] [optional: kwargs]`
* __*Example*__: `OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/Ic/60,5/V/70 dome_follow=off`

#### SKYFLAT
* __*Description*__: Morning / evening sky flats command. If `ra`, `dec` or `alt`, `az` is not given, but `seq` is
given, program will do `seq` in actual telescope position.
* __*Default*__: Coordinates epoch is J2000, dome is follow telescope, 
mirror covers are not opening/closing, dome cover is not opening/closing, **dithering is ON**
(see `dither` kwarg), tracking is ON. 
* __*Args*__: `object_name`, `ra`, `dec`
* __*Kwargs*__: `seq`, `az`, `alt`, `tracking`, `dither`, `dome_follow`, `mirror_cover`, `epoch`, `observer`, `uobi`
* __*Syntax*__: `SKYFLAT [optional: object_name] [optional: ra] [optional: dec] [optional: kwargs]`
* __*Example*__: `SKYFLAT Z022 18:58:14.75 17:21:39.29 seq=5/Ic/5,5/V/a`
!!! note
    Please notice about SKYFLAT **auto exposure option** - see `seq`.

#### DOMEFLAT
* __*Description*__: Dome flat command. The command get dome position and screen position from settings.
* __*Default*__: Mirror covers are not opening/closing, dome cover is not opening/closing, screen ligth is 
not on/off, tracking is OFF.
* __*Args*__: `object_name`
* __*Kwargs*__: `seq`, `screen_light`, `mirror_cover`, `observer`, `uobi`
* __*Syntax*__: `DOMEFLAT [optional: object_name] [optional: kwargs]`
* __*Example*__: `DOMEFLAT W1 seq=5/Ic/5,5/V/5 screen_light=auto`

#### FOCUS
* __*Description*__: Focusing command.  
If `ra`, `dec` or `alt`, `az` is not given, but `seq` it is, program will do `seq` in actual telescope position.
* __*Default*__: Coordinates epoch is `2000`, dome is follows telescope, 
mirror covers are not opening/closing, dome cover is not opening/closing, tracking is ON.
* __*Args*__: `object_name`, `ra`, `dec`
* __*Kwargs*__: `seq`,`az`, `alt`, `pos`, `auto_focus`, `tracking`, `dome_follow`, `mirror_cover`, `epoch`, `observer`, `uobi`
* __*Syntax*__: `FOCUS [optional: object_name] [optional: ra] [optional: dec] [pos] [seq] [optional: kwargs]`
* __*Example*__: `FOCUS RR1 18:58:14.75 17:21:39.29 pos=15200/100 seq=5/Ic/3 auto_focus=on`

#### DARK
* __*Description*__: This command is executing dark exposures.
* __*Default*__:  Mirror covers are not opening/closing, dome cover is not opening/closing.
* __*Args*__: `object_name`
* __*Kwargs*__: `seq`, `mirror_cover`, `observer`, `uobi`
* __*Syntax*__: `DARK [optional: object_name] [optional: kwargs]`
* __*Example*__: `DARK FF_Aql22 seq=5/Ic/60,5/V/70`

#### ZERO
* __*Description*__: This command is executing zero exposures.
* __*Default*__:  Mirror covers are not opening/closing, dome cover is not opening/closing.
* __*Args*__: `object_name`
* __*Kwargs*__: `seq`, `mirror_cover`, `observer`, `uobi`
* __*Syntax*__: `ZERO [optional: object_name] [optional: kwargs]`
* __*Example*__: `ZERO FF_Aql22 seq=5/Ic/0,5/V/0 mirror_cover=close`

#### SNAP
* __*Description*__: Command snap behaves like `OBJECT` with difference that exposure is not saved to 
data storage. Snap image is always one and is saved locally. Use this command for tests.
If `ra`, `dec` or `alt`, `az` is not given, but `seq` it is, program
will do `seq` in actual telescope position. This command is possible to dithering.
* __*Default*__:  Coordinates epoch is `2000`, dome is follow telescope, 
mirror covers are not opening/closing, dome cover is not opening/closing, dithering is `OFF`, tracking is `ON`.
* __*Args*__: `object_name`, `ra`, `dec`
* __*Kwargs*__: `seq`, `az`, `alt`, `tracking`, `dither`, `dome_follow`, `mirror_cover`, `epoch`, `observer`, `uobi`
* __*Syntax*__: `SNAP [optional: object_name] [optional: ra] [optional: dec] [optional: kwargs]`
* __*Example*__: `SNAP FF_Aql 18:58:14.75 17:21:39.29 seq=1/Ic/60 dome_follow=off`

#### STOP
* __*Description*__: This command is stops program execution.
* __*Kwargs*__: `observer`, `uobi`
* __*Syntax*__: `STOP [optional: kwargs]`
* __*Example*__: `STOP`

#### WAIT
* __*Description*__: This command is waits for time event. 
* __*Kwargs*__:  `wait`,`ut`, `wait_sunrise`, `wait_sunset`, `observer`, `uobi`
* __*Syntax*__: `WAIT [kwargs]`
* __*Example*__: `WAIT ut=22:23:33`

### Args

#### object_name
* __*Description*__: This field is identifying object and is saved to fits header. 
* __*Syntax*__: `[object name text without space]`
* __*Example*__: `FF_Aql`
!!! note
    In future will be used to get object coordinates.

#### ra
* __*Description*__: Right ascension coordinate, use only with `dec`, put in hourangle. 
* __*Syntax*__: `[right ascension in hourangle]`
* __*Range*__: `[00:00:00, 23:59:59.99]`
* __*Example*__: `18:58:14.75`
!!! warning
    Do not use with `alt`, `az` kwargs.

#### dec
* __*Description*__: Declination coordinate, use only with `ra`, put in sexagesimal.
* __*Syntax*__: `[declination in degrees sexagesimal]`
* __*Range*__: `[-90:00:00, 90:00:00]`
* __*Example*__:  `-17:21:39.29`
!!! warning
    Do not use with `alt`, `az` kwargs.

### Kwargs

#### seq
* __*Description*__: Sequence of sub-exposures with designated filters. Filter names must match with available filters in
filter wheel. Time of sub-exposure is set in seconds. **Skyflats auto exposure option:** to avoid
manual `SKYFLATS` exposure calculation, set "a" option instead exposure time in `seq` (example: `seq=5/r/a,5/V/a`).
The `SKYFLAT` command will wait to the proper light conditions to do exposures and will end if conditions will be
not satisfying anymore. How it works: exposure procedure is starting but image is not taken -> camera is making short
exposures (like 1s) and calculate mean (doing that many times if necessary) -> if mean reach set parameters 
camera making image normal way but with exposure time calculated from measure -> procedure starts again on another image. 
**Sequence multiplier:** to avoid putting same sequences in one command like `seq=1/V/1,1/r/1,1/V/1,1/r/1` there is 
possibility to use syntax `seq=2x(1/V/1,1/r/1)`.
* __*Syntax*__: `seq=[number of subexposures 1]/[filter name 1]/[time of subexposures 1],[number of subexposures 2]/[filter name 2]/[time of subexposures 2],...`
* __*Syntax OPTION*__: `seq=[number of subexposures 1]/[filter name 1]/a,...`
* __*Syntax OPTION*__: `seq=[number of repetitions]x([number of subexposures 1]/[filter name 1]/[time of subexposures 1],...)`
* __*Example*__: `seq=10/Ic/20,10/V/30,5/V/50`

#### alt
* __*Description*__: Telescope altitude, use only with `az` kwarg, put in degrees.
* __*Syntax*__: `alt=[altitude in degrees]`
* __*Range*__: `[0, 90]`
* __*Example*__: `alt=75.0`
!!! warning
    Do not use with `ra`, `dec` args.

#### az
* __*Description*__: Telescope azimuth, use only with `alt` kwarg, put in degrees.
* __*Syntax*__: `alt=[azimuth in degrees]`
* __*Range*__: `[0, 360]`
* __*Example*__: `az=182.3`
!!! warning
    Do not use with `ra`, `dec` args.

#### tracking
* __*Description*__: This kwarg is overwriting default tracking on/off.
* __*Syntax*__: `tracking=[values]`
* __*Values*__: `[on, off]`
* __*Example*__: `tracking=off`

#### dome_follow
* __*Description*__: This kwarg is overwriting default dome following telescope position. If `off` dome is not follow.
* __*Syntax*__: `dome_follow=[values]`
* __*Values*__: `[off]`
* __*Example*__: `dome_follow=off`

#### mirror_cover
* __*Description*__: This kwarg is operating mirror covers, witch from safety reason **are OFF** by default. The `auto` value
is opening covers on begin or closing on begin (for dark/zero) if necessary but not closing at end. The `open` value
opens at begin. The `close` value is closing at end and on begin (for dark/zero).
* __*Syntax*__: `mirror_cover=[values]`
* __*Values*__: `[auto, open, close]`
* __*Example*__: `mirror_cover=close`

#### screen_light
* __*Description*__: This kwarg is operating dome screen light and is OFF by default. The `auto` value
is turning on screen light on begin and turn the light off at end, `on` is turn the light on begin and `off`
is turn the light off at end. This kwarg is working only with `DOMEFLAT`.
* __*Syntax*__: `screen_light=[values]`
* __*Values*__: `[auto, on, off]`
* __*Example*__: `screen_light=on`

#### dither
* __*Description*__: This kwarg is operating dithering option. 'Dither mode' is the way how dithering works, 
for example 'basic' mode do random movement +/-'dithering distance' in both ra and dec direction. 'Dithering distance'
value is the maximum step given in minutes coordinate. 'Dithering every exposure' value is describing how often
dithering will happen (value 1 mean that dithering is working on every exposure, value 2 mean dithering is working on
first, third, fifth exposure). Syntax OPTION: `dither=off` - turning off the dithering.
* __*Syntax*__: `dither=[dither mode]/['dithering every exposure' in integer]/['dithering distance' in minutes float]`
* __*Syntax OPTION*__: `dither=off`
* __*Values*__: `dither mode = [basic]`
* __*Example*__: `dither=basic/2/1`
!!! note
    In future there will be more 'Dither mode' added.

#### epoch
* __*Description*__: This kwarg is operating `ra`, `dec` coordinate epoch. 
* __*Syntax*__: `epoch=[epoch]`
* __*Values*__: `epoch = [..., 1975, 2000, 2025, 2050, ...]`
* __*Example*__: `epoch=2025`

#### pos
* __*Description*__: This kwarg is operating focusing procedure. 'Focusing target' is approximately focusing target.
'Focusing step' is step value used to create focusing curve. For example `pos=15700/100` value will create focuser steps
`[15500, 15600, 15700, 15800, 15900]`.
* __*Syntax*__: `pos=['focusing target' in integer]/['focusing step' in integer]`
* __*Example*__: `pos=15700/100`

#### auto_focus
* __*Description*__: This kwarg is operating focusing. With `pos` kwarg is composing automatically focusing procedure.
Sentence auto_focus=off is equal to not mention this kwarg at all.
* __*Syntax*__: `auto_focus=[on, off]`
* __*Example*__: `auto_focus=on`

#### wait
* __*Description*__: This kwarg is specifying how many second program will wait.
* __*Syntax*__: `wait=[seconds float]`
* __*Example*__: `wait=320.0`

#### ut
* __*Description*__: This kwarg is specifying UT time to wait.
* __*Syntax*__: `ut=[UT time]`
* __*Example*__: `ut=0:03:23`

#### wait_sunrise
* __*Description*__: This kwarg is specifying sunrise position to wait.
* __*Syntax*__: `wait_sunrise=[float sun alt position]`
* __*Example*__: `wait_sunrise=-8`

#### wait_sunset
* __*Description*__: This kwarg is specifying sunset position to wait.
* __*Syntax*__: `wait_sunset=[float sun alt position]`
* __*Example*__: `wait_sunset=-18`

#### uobi
* __*Description*__: This kwarg is identification number used by software.
* __*Syntax*__: `uobi=[integer number]`
* __*Example*__: `uobi=12222351`

#### observer
* __*Description*__: This kwarg is specifying observers name that will be added to fits header. 
* __*Syntax*__: `observer=[observers name]`
* __*Example*__: `observer=JonArt_MikeSmith_AnnR`
!!! note
    Please don't use space in observers name.

### Comments
We will use "line comments"  indicated by the hash `#` sign. Parser will ignore `#` and anything which follows it to the
end of the line. 
!!! note
    Please use space after # (example: # OBJECT seq=3/r/5).
