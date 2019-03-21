
			MOOSH

GENERAL OVERVIEW

The aim of Moosh is to provide a complete set of tools to compute all
the more or less usual optical properties of any multilayered
structure: reflexion, transmission, absorption spectra, as well as
gaussian beam propagation or guided modes. It can be seen as a
semi-analytic (making it light and fast) solver for Maxwell's
equations in multilayers. It is written in Octave/Matlab, available on
Github and based on scattering matrices, making it perfectly
stable. This software is meant to be extremely easy to (re)use, and
could prove useful in many research areas like photovoltaics,
plasmonics and nanophotonics, as well as for educational purposes for
the high number of physical phenomenon it can illustrate.

The implementation of Moosh has actually much to do with a real army
knife. A geometrical and an electromagnetic descriptions of the
structure are made in a central file structure.m in such a way that it
is very simple to describe even very complex multilayers. Many
programs are then available that are able to compute various optical
properties of the structure described in structure.m.

DESCRIBING THE STRUCTURE

In file "structure.m" the different media that can be considered are
first defined. Each of them is characterized by a permittivity and a
permeability, that are contained in two vectors called epsilon and mu,
respectively.

Example : epsilonn=[1,2.25]; muu=[1,1]; Here, epsilon and mu are
lists. The first environment has a permittivity epsilon(1)=1 and a
permeability mu(1)=1. The second one has a permittivity
epsilon(2)=2.25 and a permeability mu(2)=1. So, the first one
corresponds to air here and the second one to glass with an index of
1.5.

The structure is then described by specifying how the media are
stacked. Light is assumed to come from above, and the media are given
in the very order light will encounter them.

Example: typee=[1,2,1,2,1];

Here, there are 5 layers : the first is constituted by air, the second
by glass, the third by air, the forth again by glass and the last one
by air.

The height of each layer has then to be specified.

Example : hauteur=[20*600,600,600,600,20*600]; The first layer has a
height of 20*600, the second a height of 600, the third 600, the forth
600 and the fifth 20*600. You must take care to list the heights in
the same order as in typee.

The polarization has to be chosen there too, even if some functions
can override this.  Example : pol=1; In this case, the polarization
will be H// (TM). If pol=0, the polarization is E// (TE).

MAIN PROGRAMS

Each of these programs requires the user to specify a few parameters
concerning either the physical situation, or the window that must be
represented. These parameters are all located at the beginning of each
program, except for the vizualisation parameters, that can be found at
the end.

* Reflection, transmission and absorption

Angular.m : To get the reflection, transmission and absorption
coefficients of the structure described in "structure.m" as a function
of the angle.

Spectrum.m : To get the reflection, transmission and absorption
coefficients of the structure described in "structure.m" as a function
of the wavelength.

* Beam propagation

Beam.m : To get a picture of the intensity of the field (electric or
magnetic field, depending on the polarization) when the structure is
illuminated with a gaussian beam characterized by a wavelength, a
waist (typical width) and an incidence angle.

* Photovoltaics

Photo.m : Computes the theoretical short-circuit current for a given
structure containing an active layer, as well as the conversion
efficiency.  This program uses an AM 1.5 spectrum (am1_5.m in data/),
as well as absorption.m to compute the absorption.

* Guided modes

Guidedmodes.m : This code search transverse guided modes of the
structure to build a vector containing all propagation constants of
every guided modes of the structure. The resulting propagation
constantes are stored in a "modes" variable.

Profile.m : Profile must be called after Guidedmodes.m, typically with
Profile(modes(n),lambda) to compute the profile of the n-th mode that
has been found by Guidedmodes.m.

Map.m : The "dispersion.m" function allows to see the solutions of 
the dispersion relation of a multilayered structure in the complex
plane as zeros of "dispersion.m". Any zero of this function 
corresponds to a guided mode. 

* Punctual sources

Green.m : Computes the Green function when a punctual source is placed
in the structure.

HOW TO RUN BEAM.M

Beam.m is probably the most useful tool to visualize a lot of
phenomenon, so that its use is described here.

A lot has to be specified about the incoming beam: the incidence
angle, the wavelength, the waist (defined here as the horizontal waist
at the very top of the structure).  The rest of the parameters deals
with the simulation window: its size, the number of points for which
the field is computed (nx pixels in the horizontal direction, ny in
the vertical one), the position of the beam inside the simulation
window (C=0 means to the left, C=1 means to the right, and C=0.5 is in
the very middle of the simulation window).

Example : lambda=600; d=100*lambda; w=10*lambda; theta=56.55*pi/180;
	   C=0.2; nx=d/50 ny=floor(hauteur/50)
       
Here, the wavelength is 600 nm, the incident angle 56.55°, the size of
the simulation window 60000 nm, the waist of the beam 6000 nm, the
incident beam is placed at the top at 20% of the simulation window
from the left corner.

LICENSE

This program is distributed under the GNU General Public License 2.0

CONTACT

For any question, please contact antoine.moreau@univ-bpclermont.fr.
