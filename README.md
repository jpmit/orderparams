orderparams
===========

James Mithen  
j.mithen@surrey.ac.uk

orderparams is a C++ code for computing structural properties (or
'order parameters') of systems of particles.  It is specifically
designed to analyse structural properties of configurations that
result from simulations of crystallisation.  For example, it will
attempt to quantify which particles are in a crystalline environment
using methods that have been designed specifically for this purpose.

There are two executables that can be built.  The main one is
'orderparams', which computes a number of structural properties (or
'order parameters') of a system of particles, where the positions of
the particles are given in the XYZ file format.

The other executable, 'ldtool', again takes the positions of the
particles in the XYZ file format as an input and uses the
Lechner-Dellago method for identifying crystalline particles.

See USAGE below for details of how to run the executables.

REQUIREMENTS
-------------

* C++ compiler (tested with gcc 4.6.3)
* Boost C++ libraries (tested with version 1.46)
* GNU scientific libraries (GSL)

COMPILING
-----------

The Makefile has two possible targets (i.e. two executables can be
built).  The default one is the main executable, 'orderparams'.  To
compile 'orderparams', navigate to the root directory and type

    $ make

To compile the other executable 'ldtool', navigate to the root
directory and type

    $ make ldtool

USAGE
--------

An example input file, 'params.out' is contained in the examples/
directory.  To run both executables for the example provided, navigate
to the root directory and type

    $ ./orderparams examples/params.out

and 
    
    $ ./ldtool examples/params.out


OUTPUT OF orderparams
------

orderparams outputs the following order parameters (42 in total) for
the XYZ configuration specified in the input file (in the examples
directory, this is the file named 'params.out'.  The 'Name' field
below is how the output appears in the terminal, and the 'Description'
field gives a brief description:

<table>
  <tr>
    <th>Number</th><th>Name</th><th>Description</th>
  </tr>
  <tr>
    <td>1</td><td>N_ld</td><td>Size of largest crystalline cluster using LD method</td>
  </tr>
  <tr>
    <td>2</td><td>N_tf</td><td>Size of largest crystalline cluster using TF method</td>
  </tr>
  <tr>
    <td>3</td><td>n_bccLD</td><td>Fraction of bcc pars in LD cluster</td>
  </tr>
  <tr>
    <td>4</td><td>n_bccTF</td><td>Fraction of bcc pars in TF cluster</td>
  </tr>
  <tr>
    <td>5</td><td>n_fccLD</td><td>Fraction of fcc pars in LD cluster</td>
  </tr>
  <tr>
    <td>6</td><td>n_fccTF</td><td>Fraction of fcc pars in TF cluster</td>
  </tr>
  <tr>
    <td>7</td><td>n_hcpLD</td><td>Fraction of hcp pars in LD cluster</td>
  </tr>
  <tr>
    <td>8</td><td>n_hcpTF</td><td>Fraction of hcp pars in TF cluster</td>
  </tr>
  <tr>
    <td>9</td><td>n_icosLD</td><td>Fraction of icosahedral pars in LD cluster</td>
  </tr>
  <tr>
    <td>10</td><td>n_icosTF</td><td>Fraction of icosahedral pars in TF cluster</td>
  </tr>
  <tr>
    <td>11</td><td>Q6clusLD</td><td>Average Q6 of LD cluster</td>
  </tr>
  <tr>
    <td>12</td><td>Q6clusTF</td><td>Average Q6 of TF cluster</td>
  </tr>
  <tr>
    <td>13</td><td>Q4clusLD</td><td>Average Q4 of LD cluster</td>
  </tr>
  <tr>
    <td>14</td><td>Q4clusTF</td><td>Average Q4 of TF cluster</td>
  </tr>
  <tr>
    <td>15</td><td>N_sLD</td><td>Number of liquid like particles with at least one neighbour in LD cluster</td>
  </tr>
  <tr>
    <td>16</td><td>N_sTF</td><td>Number of liquid like particles with at least one neighbour in TF cluster</td>
  </tr>
  <tr>
    <td>17</td><td>N_lLD</td><td>Total number of crystalline connections (or links) for all liquid-like particles with at least one neighbour in LD cluster</td>
  </tr>
  <tr>
    <td>18</td><td>N_lTF</td><td>Total number of crystalline connections (or links) for all liquid-like particles with at least one neighbour in TF cluster</td>
  </tr>
  <tr>
    <td>19</td><td>Q6N_sLD</td><td>Average q6 of liquid-like particles with at least one neighbour in LD cluster</td>
  </tr>
  <tr>
    <td>20</td><td>Q6N_sTF</td><td>Average q6 of liquid-like particles with at least one neighbour in TF cluster</td>
  </tr>
  <tr>
    <td>21</td><td>Q4N_sLD</td><td>Average q4 of liquid-like particles with at least one neighbour in LD cluster</td>
  </tr>
  <tr>
    <td>22</td><td>Q4N_sTF</td><td>Average q4 of liquid-like particles with at least one neighbour in TF cluster</td>
  </tr>
  <tr>
    <td>23</td><td>Rbar_g,1LD</td><td>Smallest eigenvalue of gyration tensor of LD cluster</td>
  </tr>
  <tr>
    <td>24</td><td>Rbar_g,1TF</td><td>Smallest eigenvalue of gyration tensor of TF cluster</td>
  </tr>
  <tr>
    <td>25</td><td>Rbar_g,2LD</td><td>Middle eigenvalue of gyration tensor of LD cluster</td>
  </tr>
  <tr>
    <td>26</td><td>Rbar_g,2TF</td><td>Middle eigenvalue of gyration tensor of TF cluster</td>
  </tr>
  <tr>
    <td>27</td><td>Rbar_g,3LD</td><td>Largest eigenvalue of gyration tensor of LD cluster</td>
  </tr>
  <tr>
    <td>28</td><td>Rbar_g,3TF</td><td>Largest eigenvalue of gyration tensor of TF cluster</td>
  </tr>
  <tr>
    <td>29</td><td>Rbar_gLD</td><td>Square of 'radius of gyration' for LD cluster</td>
  </tr>
  <tr>
    <td>30</td><td>Rbar_gTF</td><td>Square of 'radius of gyration' for LD cluster</td>
  </tr>
  <tr>
    <td>31</td><td>R_g,zLD</td><td>(3,3) element of non-diagonalized gyration tensor of LD cluster</td>
  </tr>
  <tr>
    <td>32</td><td>R_g,zTF</td><td>(3,3) element of non-diagonalized gyration tensor of TF cluster</td>
  </tr>
  <tr>
    <td>33</td><td>R_g,1LD</td><td>Smallest eigenvalue of top-diagonalised gyration tensor of LD cluster</td>
  </tr>
  <tr>
    <td>34</td><td>R_g,1TF</td><td>Smallest eigenvalue of top-diagonalised gyration tensor of TF cluster</td>
  </tr>
  <tr>
    <td>35</td><td>R_g,2LD</td><td>Middle eigenvalue of top-diagonalised gyration tensor for LD cluster</td>
  </tr>
  <tr>
    <td>36</td><td>R_g,2TF</td><td>Middle eigenvalue of top-diagonalised gyration tensor for TF cluster</td>
  </tr>
  <tr>
    <td>37</td><td>s_bcc</td><td>Fraction of bcc particles in entire system (excluding surface particles)</td>
  </tr>
  <tr>
    <td>38</td><td>s_fcc</td><td>Fraction of fcc particles in entire system (excluding surface particles)</td>
  </tr>
  <tr>
    <td>39</td><td>s_hcp</td><td>Fraction of hcp particles in entire system (excluding surface particles)</td>
  </tr>
  <tr>
    <td>40</td><td>s_icos</td><td>Fraction of icosahedral particles in entire system (excluding surface particles)</td>
  </tr>
  <tr>
    <td>41</td><td>Q6</td><td>Average q6 of all particles in system (excluding surface particles)</td>
  </tr>
  <tr>
    <td>42</td><td>Q4</td><td>Average q4 of all particles in system (excluding surface particles)</td>
  </tr>
</table>


OUTPUT OF ldtool
------

ldtool will output an XYZ file (to stdout) with the same particle
positions as specified in the input file (in the examples directory,
the input file is params.out), and with symbols that identify the
classification 'type' of the particle.  This uses the Lechner-Dellago
approach to identifying crystalline particles.  The symbols map to
particle types as given below.  Also given is the colour that the
particle is depicted as when viewed using the visualisation software
Jmol.

<table>
  <tr>
    <th>Symbol</th><th>Type</th><th>Jmol colour</th>
  </tr>
  <tr>
    <td>S</td><td>Fcc</td><td>Yellow</td>
  </tr>
  <tr>
    <td>P</td><td>Hcp</td><td>Orange</td>
  </tr>
  <tr>
    <td>F</td><td>Bcc</td><td>Green</td>
  </tr>
  <tr>
    <td>O</td><td>Surface</td><td>Red</td>
  </tr>
  <tr>
    <td>N</td><td>Liquid</td><td>Blue</td>
  </tr>
  <tr>
    <td>O</td><td>Icosahedral</td><td>Pink</td>
  </tr>
</table>

TODO
-----------

* improve docs (this file and source files)

* make error handling more comprehensive

* remove dependency on GSL library (use Boost instead?)

* simplify Makefile

* improve reading parameters from params.out (?). Could possibly read
  from params.pkl with interface to python code.

* in box.h make box lengths private and create method(s) for changing
  size

* add more order parameters?
