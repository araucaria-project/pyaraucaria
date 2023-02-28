This file was copied form https://github.com/araucaria-project/doc-oca-tic-toi and extended.

# Observations Plan Files for OCA observatory

At the current stage, all solutions and architectural decisions presented here have the status of preliminary proposals 
and are intended as an invitation to the discussion.

## Objectives

After evaluation of pros and cons, we decided to introduce new file format for observation plans, which will be used by the new (2022+) OCA observatory software.
As, in some simplification, observations plan is a sequence of command, the format of Observations Plan file is basically a syntax of kind of programming language.
Because we want the new format to be human-readable, -editable and -manageable, we decide to introduce new language instead of using e.g. JSON syntax. 
New format should be somehow similar to ols obsplan files for the **IRIS** and **V16** telescopes and in general to be readable and understandable by astronomers without 
need to RTFM (*read the fantastic manual*).

We want to keep the syntax well-defined and intuitively but uniquely transformable to python dictionary object (and therefore to JSON, YAML etc.).
Moreover, we want to provide any developers with python package for parsing anf formatting those files.

## The Syntax

We will introduce the proposed syntax step by step showing how to fulfill collection of specific demands

### Simplicity - Single Command
We want the new language to act as internal language of the software, so the simplest form is just single command to be executed rather than full-fledged observation plan.
E.g. we want following text to be proper "program":
```
  WAIT t=20
```
or, a bit more complicated single command
```
  OBJECT HD193901 20:23:35.8 -21:22:14.0 seq = 5/I/60,5/V/70
```
Both of single command lines are understandable for OCA astronomers 
(wait for 20 seconds, make photos of object HD193901 and ra/dec coordinates 20:23:35.8 -21:22:14.0 five times for 50 seconds with filter *I* 
and five times for 70 seconds with filter *V*). Syntactical break down of the last command goes as follows:
* command: `OBJECT`
* positional arguments: `HD193901`, `20:23:35.8`, `-21:22:14.0` - some from the tail may be optional
* keyword arguments: argument named `sequence` of value `5/I/60,5/V/70`. all keyword arguments can be optional

One may say, that it would be more familiar to write just `OBJECT HD193901 20:23:35.8 -21:22:14.0 5/I/60,5/V/70`, indeed 
it's simpler, but we propose coordinates arguments to be optional (we have database of filed position) so one could also write just
```
  OBJECT HD193901 sequence = 5/I/60,5/V/70
```
but we can consider another semantics of the `OBJECT` command.


## Sequences
Sequence is a number of commands each written in unique, single line e.g:
```
WAIT ut=16:00
ZERO seq=15/I/0
DARK seq=10/V/300,10/I/200
DOMEFLAT seq=7/V/20,7/I/20
DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
WAIT sunset=-12
SKYFLAT alt=60:00:00 az=270:00:00  seq=10/I/20,10/V/30 
SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
WAIT t=600
FOCUS NG31 12:12:12 20:20:20
OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70
OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20 focus=+30
```

Sometimes, more often than we would like, we have to ask telescope to (re)start sequence from specific instruction 
instead of from the beginning. For such situation, we need some method to refer to this instruction. We see two simplest methods:
we can refer to particular line number, or we can introduce explicit labels. As we want to keep labels optional,
we will use semicolon at the end of the label to avoid confusion label with command. 

The previous example rewritten with sample labels for each line:
```
START: WAIT ut=16:00
00100: ZERO seq=15/I/0
00110: DARK seq=10/V/300,10/I/200
00120: DOMEFLAT seq=7/V/20,7/I/20
00130: DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
SUNSET: WAIT sunset=-12
00150: SKYFLAT alt=60:00:00 az=270:00:00  seq=10/I/20,10/V/30 
00160: SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
00170: WAIT t=600
00100: FOCUS NG31 12:12:12 20:20:20
OB01:  OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
OB02:  OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70
OB03:  OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20 focus=+30
```

Introduction of labels brings as to something like below as a syntax of single command:
In general the syntax of single command line should be like that:
```
    [LABEL:] <COMMAND> [<POSITIONAL_ARG1> [<POSITIONAL_ARG1>...]] [KW_ARG1_NAME=KW_ARG1_VAL [KW_ARG2_NAME=KW_ARG2_VAL]...]
```

### Syntax Details
#### Comments
We will use "line comments"  indicated by the hash `#` sign. Parser will ignore `#` and anything which follows it to the
end of the line.

#### Typing
Parser will try to derive type of any parameter:
* First parser will try to cast any parameter to `int`,
* If it fails, then parser will try to cast parameter to `float`,
* If it fails, then parser treats parameter as `str`.

Because space is an argument separator, string parameters can be put in optional quotation marks `"` or apostrophes `'`
(which will not became parts of parameter value). One can also use python escape characters e.g. `\"` for quotation 
mark inside string or `\n` for new line.

### Possible Extensions
In future variables and arithmetical expression can be introduced, so our suggestion is to limit allowed characters
in unquoted string values (e.g. forbid arithmetical operators and '$' sign which will be used to denote variable name).

## Non-sequential Execution
There are situations, when sequential execution of an Observation Program is not enough. Below we discuss such 
circumstances.
### Circumstances Breaking Sequential Execution
#### 1. Restarting script from some point
After break of script execution, caused by harsh weather, failure or another reason, we are restarting script from
specific line rather than from the beginning. However, there may be some instructions which have to be executed on
each start or restart, e.g. bringing camera to operational temperature (it may not be real-life example, but good for 
illustration purpose). Therefore, we need some section with commands executed always on restarting.

#### 2. Scheduled Observation
We may have observations to be performed exactly at specific time (e.g. during eclipse od EBS). We need therefore
a section (or at least command) which will be scheduled for that ime and will break the sequence when the time comes.

#### 3. Periodic Calibrations
We may need e.g. to refocus camera periodically, with constant delay, or, even worst, on specific change of the
environment temperature. Similarly to "Scheduled Observation", it breaks sequence, and may be hard to determine
the exact moment of this interruption.

#### 4. End of Night Sequence
Commands for morning flat-fielding nad closing telescopes used to be added at the end of sequence. It worked well,
if the execution time of whole program was well established. Sometimes we were losing some time waiting for the dawn,
with conservative observation plan leaving safe gap after last observation. But having scheduling mechanisms anyway,
we can imagine that the dawn sequence is just scheduled for specific moment (expressed as time or sun's altitude).

### Resolution
For the *executor* (software executing Observation Plan), it's possible to schedule non-sequential parts of program.
Anyway we have to introduce some syntax for those parts. It looks like, it's unavoidable to introduce som kind of
blocks of instructions. Blocks of instructions, are of course present in any modern programming languages.

Here are some requests for the syntax of blocks of instructions for our "language":
* We are trying to avoid use of brackets od any kind, so it's better to do not use any brackets for indication of
the boundaries ot the block (no C++, no Java).
* It would be better to do not relay on white characters (no Python). It should be however possible to use indentation
for visually mark the blocks, but this should be not obligatory.
* Begin and end of the block instruction should fit in one line. (e.g. no "`BLOCK` `BEGIN`" but rather 
single "`BEGINBLOCK`")  

We would suggest something like this
```
BEGINSEQUENCE execute_at_time=16:00
    ZERO seq=15/I/0
    DARK seq=10/V/300,10/I/200
    DOMEFLAT seq=7/V/20,7/I/20
    DOMEFLAT seq=10/str_u/100 domeflat_lamp=0.7
ENDSEQUENCE

BEGINSEQUENCE execute_at_time=02:21:43 priority=+30  # scheduled obs
    OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=2/I/60,2/V/70
ENDSEQUENCE

BEGINSEQUENCE execute_periodically=02:00 priority=+10
    FOCUS NG31 12:12:12 20:20:20
ENDSEQUENCE

BEGINSEQUENCE execute_at_dusk=-12
    SKYFLAT alt=60:00:00 az=270:00:00  seq=10/I/20,10/V/30 
    SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
    WAIT t=600
    FOCUS NG31 12:12:12 20:20:20
    OBJECT HD193901 20:23:35.8 -21:22:14.0 seq=1/V/300
    OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70
    OBJECT V496_Aql 19:08:20.77 -07:26:15.89 seq=1/V/20 focus=+30
    # etc ...
ENDSEQUENCE

BEGINSEQUENCE execute_at_dawn=-6 priority=+10
    SKYFLAT alt=60:00:00 az=270:00:00  seq=10/I/20,10/V/30 
    SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
    SKYFLAT alt=60:00:00 az=270:00:00  seq=10/I/20,10/V/30 
    SKYFLAT seq=10/I/0,10/V/0 skyflat_adu=30
ENDSEQUENCE

BEGINSEQUENCE execute_at_dawn=+2 priority=+100
    PARK 
    DOMECLOSE
ENDSEQUENCE
```
Please note, that morning flats will break execution of observations (`prority=+10`) and telescope closing
will similarly interrupt dawn flat-fielding (`prority=+100`). 

## Parser
We will implement Observation Plan parser (and maybe also formatter) as independent module, probably in 
`https://github.com/araucaria-project/pyaraucaria` repository. It will allow any party to handle this file format
in common way.

we suggest, that the implementation of the parser will base on the 
[lark python package](https://github.com/lark-parser/lark) which is the most known and used abstract parser. The
advantage of using such parser is also fact that we will have to have formal EBNF syntax definition.

The parser will translate Observation Plan to python `dict` object. E.g. the last example should produce following 
dictionary:
```python
{'commands': [
    {
        'command': 'SEQUENCE',
        'kwargs': {'execute_at_time': '16:00'},
        'commands': 
        [
            {
                'command': 'ZERO',
                'kwargs': {'seq': '15/I/0'}
            },
            {
                'command': 'DARK',
                'kwargs': {'seq': '0/V/300,10/I/200'}
            },
            {
                'command': 'DOMEFLAT', 
                'kwargs': {'seq': '7/V/20,7/I/20'}
            },
            {
                'command': 'DOMEFLAT',
                'kwargs': {'seq': '10/str_u/100','domeflat_lamp': 0.7}
            }
        ]  
    },
    {
        'command': 'SEQUENCE',
        'kwargs': {'execute_at_time': '02:21:43','priority': 30},
        'commands': [
            {
                'command': 'OBJECT',
                'args': ['FF_Aql', '18:58:14.75', '17:21:39.29'],
                'kwargs': {'seq': '2/I/60,2/V/70'}
            }
        ]
    }
]}
```
etc...

Beside `command`, `commands`, `args` and `kwargs`, we consider following tags for command:
* `ln` - line number, for debugging and restarting,
* `label` if provided
* `message` may be added on latter phases of processing - what to display on execution monitor
and others...

All those names are still subject to change (e.g. `command` and `commands` can be confusing).

# We're getting close to a solution
This section summarizes all of the above considerations for a final solution.
## Main commands
Here are main commands used in observation plans:
#### OBJECT
* Description: This command will unpark mount, turn on tracking, slew to sky object co-ordinates. Also, will check if 
dome is open, if not it will open and slew dome to az telescope position and finally sync moves. OBJECT command strongly goes
together with "seq" kwarg.
* Args: object name - this name will be used to put into fits header also in future can be used to get object
co-ordinates from catalogue; ra, dec - object co-ordinates, in future optional.
* Kwargs: seq, focus, az, alt
* Syntax: OBJECT [object name] [optional: ra] [optional: dec] [optional: kwargs]
* Example: OBJECT FF_Aql 18:58:14.75 17:21:39.29 seq=5/I/60,5/V/70

#### ZERO
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### DARK
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### DOMEFLAT
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### SKYFLAT
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### WAIT
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### FOCUS
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### BEGINSEQUENCE
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 
#### ENDSEQUENCE
* Description: 
* Args:
* Kwargs:
* Syntax:
* Example: 

## Kwargs
#### seq
* Description: Sequence of sub-exposures with designated filters. Filter names must match with available filters in
filter wheel. Time of sub-exposure is set in seconds. If main command will prompt with "seq", camera immediately will 
set to parameters and filter wheel will rotate to proper position.
* Syntax: seq=[number of subexposures 1]/[filter name 1]/[time of subexposures 1], 
[number of subexposures 2]/[filter name 2]/[time of subexposures 2],...
* Example: seq=10/I/20,10/V/30,5/V/50
#### focus
* Description:
* Syntax:
* Example:
#### skyflat_adu
* Description:
* Syntax:
* Example:
#### priority
* Description:
* Syntax:
* Example:
#### execute_at_dawn
* Description:
* Syntax:
* Example:
#### execute_at_time
* Description:
* Syntax:
* Example:
#### t
* Description:
* Syntax:
* Example:
#### alt
* Description:
* Syntax:
* Example:
#### az
* Description:
* Syntax:
* Example:
#### execute_periodically
* Description:
* Syntax:
* Example:
#### domeflat_lamp
* Description:
* Syntax:
* Example:
#### ut
* Description:
* Syntax:
* Example:


## References
* Lark parser package: https://lark-parser.readthedocs.io/en/latest/index.html 
* Lark Web IDE: https://www.lark-parser.org/ide/

* Regexp in python - module `re`
* Regexp Web IDE: https://regex101.com 
