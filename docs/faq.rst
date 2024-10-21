.. _faq:

Frequently asked questions
==========================


Here are some answered questions that may come to mind while reading this guide.  Many of them have been asked by real users in the past!

.. 
    how to write a faq:
    https://github.com/readthedocs/readthedocs.org/blob/11.1.2/docs/user/faq.rst?plain=1
    https://github.com/sphinx-doc/sphinx/blob/master/doc/faq.rst?plain=1

.. contents::
    :local:


Terminology dictionary
***********************

..
    https://tex.stackexchange.com/questions/121198/how-to-adjust-a-table-to-fit-on-page

.. 
    raw:: latex
    \begin{adjustbox}{width=1.2\textwidth,center=\textwidth}

.. list-table:: Different names for different simulators
    :class: longtable
    :widths: auto
    :header-rows: 1
    
    *   - Term
        - EDEN
        - `NEURON <https://neuron.yale.edu>`__
        - `GENESIS <http://genesis-sim.org/>`__
        - `Arbor <https://arbor-sim.org>`__
        - `BRIAN <https://brian2.readthedocs.io/>`__
        - Comments
    *   - Point neuron
        - ``<ComponentType extends="baseCell">``
        - ``ARTIFICIAL_CELL``
        - 
        - What's not a *cable cell*
        - What's not a ``SpatialNeuron``
        - A neuron model with the morphology largely simplified out.
    *   - Multi-compartment cell
        - ``<cell>``
        - what's not an ``ARTIFICIAL_CELL``
        - All of them really
        - cable cell
        - ``SpatialNeuron``
        - A neuron model of multiple finely detailed *finite elements*.
    *   - finite element
        - compartment
        - compartment, segment
        - compartment
        - control volume (CV)
        - compartment
        - A piece of neuron whose state is assumed to be uniform.   
    *   - Tracing point
        - A ``<segment>``\'s ``<proximal>`` and ``<distal>``
        - ``pt3d``
        - Same as compartment
        - ``point`` or  ``segment`` for a pair
        - Same as compartment
        - A ``(x,y,z,d)`` data point for a traced neurite.
    *   - Neurite span
        - ``sao864921383``
        - section
        - 
        - branch
        - section
        - A *distinct* unbranched section of neuritic cable.
    *   - Arbitrary dependence
        - ``<VariableReference>``
        - ``POINTER``
        - ``addmsg``
        - 
        - ``linked``
        - An continuous-time dependency of a variable on *any* other variable. 
    *   - Aggregate flux
        - `Multiflux <extension_multiflux.ipynb>`__
        - ``WRITE``
        - ``addmsg``
        - ``WRITE``
        - ``summed``
        - *Distinct* flows of chemical species, supplied by various mechanisms.
    *   - Spike train
        - :ref:`\<spikeArray\> <neuroml:schema:spikearray>`, `\<EdenEventSetReader\> <extension_io.ipynb#Event-series-with-EdenEventSetReader>`__
        - `VecStim <https://www.neuron.yale.edu/phpBB/viewtopic.php?f=28&t=2117>`__
        - TBD
        - Spiking cell
        - ``SpikeGeneratorGroup`` and others
        - A disembodied source of spikes
    *   - Time series
        - `\<EdenTimeSeriesReader\> <extension_io.ipynb#Time-series-with-EdenTimeSeriesReader>`__ 
        - :any:`Vector.play() <neuron:/guide/bio_faq.rst#how-do-i-simulate-a-current-clamp-with-non-pulse-behavior>`
        - TBD
        - 
        - ``value_name(x,t)``
        - A disembodied evolving quantity

.. 
    raw:: latex
    \end{adjustbox}

.. 
    When I type a new value into a numeric field, it doesn't seem to have any effect.
    \... consider... NEXT

Using the simulator
*******************

What units are simulation results recorded in?
----------------------------------------------

The unspoken NeuroML convention is to use `SI-derived units <https://en.wikipedia.org/wiki/Fundamental_unit>`_ that are products of the seven base SI units.  Keep in mind that that the typical real-life unit for concentration is :math:`mol/L` whereas the SI unit is :math:`mol/m^3`!

If you prefer quantities to be recorded in specific units, refer to :doc:`extension_io`.


The simulation output is too large!  How can I record values every millisecond or so? 
-------------------------------------------------------------------------------------

Using `<𝙴𝚍𝚎𝚗𝙾𝚞𝚝𝚙𝚞𝚝𝙵𝚒𝚕𝚎> <extension_io.ipynb>`_.

The neuron potential recording is too large!  How can I record only when spikes happened?
-----------------------------------------------------------------------------------------

Using :ref:`<𝙴𝚟𝚎𝚗𝚝𝙾𝚞𝚝𝚙𝚞𝚝𝙵𝚒𝚕𝚎> <neuroml:userdocs:quantitiesandrecording:events>` or `<𝙴𝚍𝚎𝚗𝙴𝚟𝚎𝚗𝚝𝙾𝚞𝚝𝚙𝚞𝚝𝙵𝚒𝚕𝚎> <extension_io.ipynb>`_.

How can I have the same result, every time I run my randomised simulation?
--------------------------------------------------------------------------

EDEN does have a provision for that; set the ``seed`` attribute in your `<𝙴𝚍𝚎𝚗𝙾𝚞𝚝𝚙𝚞𝚝𝙵𝚒𝚕𝚎> <intro_neuroml.ipynb#Setting-the-randomisation-seed>`__ tag. The generated random numbers should be the same at least per model contents, and program version.

But this is *not* the end of the story!  There are ways for a small perturbation to slip through and be amplified (in principle) indefinitely by chaotic systems, like neural networks tend to be.

The ways that perturbations can slip in the numbers are allowed by the finite accuracy in out computers and the resulting round-off that takes place during calculations, and caused by carrying out the calculations differently at times.

These are some factors that can lead to calculations being run differently:

- Different (versions of) the computers' operating systems and `code generators <https://en.wikipedia.org/wiki/Compiler_(computing)>`_;
- The need for speed that pushes us to do the processing in whichever way it is convenient on each machine;
- A different order of summation, that happens they are split over multiple cores in a computer;
- The distribution of model parts between computers in a multi-machine simulation.

To study the issue and potential causes further, refer to a web search on "floating point determinism" (for example, ""`Determinism and Reproducibility
in Large-Scale HPC Systems <https://wodet.cs.washington.edu/wp-content/uploads/2013/03/wodet2013-final12.pdf>`__").

In conclusion, measures to attempt deterministic simulation may get the output to be close enough in between simulation runs, but this is regrettably not guaranteed within practical limits.

The good news is that a neural network's function should not rely on one specific pick of random samples, otherwise the random numbers wouldn't be random (being part of the model themselves).  It follows that if a model's correct behaviour appears only for a small set of random ``seed``\ s, then that is the model's fault.

If some random number sequences are of critical importance (to replicate findings exactly), consider writing down the exact waveforms and then playing them back with EDEN's `<𝙴𝚍𝚎𝚗𝚃𝚒𝚖𝚎𝚂𝚎𝚛𝚒𝚎𝚜𝚁𝚎𝚊𝚍𝚎𝚛> extension <extension_io.ipynb#Time-series-with-EdenTimeSeriesReader>`__ to NeuroML, in place of the places where ``random()`` would be used.


How should I use EDEN on a (MPI) computer cluster?
--------------------------------------------------

Just use the same `command line <https://eden-simulator.org/repo#from-the-command-line>`_ as usual on all processes launched, targeting the NeuroML file to run. Refer to `Building for MPI <https://eden-simulator.org/repo#building-for-mpi>`__ for how to build from source on a cluster.

..
    LATER Running on a MPI cluster.


NeuroML resources
*****************


Where can I find NeuroML models?
--------------------------------

Most published NeuroML models can be found in the following places:

* `NeuroML-DB <http://neuroml-db.org>`__ for a systematised database of neuron and mechanism models, and a few network models;
* `Open Source Brain <https://v1.opensourcebrain.org/projects>`__, for more diverse networks and component models in NeuroML, and also other technologies;
* NeuroML models in the `ModelDB <https://modeldb.science/modellist/154351?all_simu=true>`__.   Browse also the entire ModelDB for lots of models that could be ported to NeuroML.

See also the :doc:`NeuroML guide <neuroml:Userdocs/FindingNeuroMLModels>` on the subject.

Which other tools can I use NeuroML with?
-----------------------------------------

The official NeuroML guide provides an up-to-date list of these tools:

* :doc:`The official NeuroML tooling <neuroml:Userdocs/Software/Software>`
* :doc:`Third-party tools that can work with NeuroML <neuroml:Userdocs/Software/SupportingTools>`



Setting up and running models
*****************************

What is "discretisation" or "compartments" and how does it affect my model?
---------------------------------------------------------------------------

Refer to the `"Discretisation" section <intro_spatial.ipynb#Simulation-aspect:-Discretisation-into-compartments>`__ and try the `'simple cable' exercise <intro_spatial.ipynb#Example:-Modelling-and-simulating-a-stretch-of-neural-cable>`__ exercise, on the chapter on spatially-modelled neurons.  What usually happens is:

* Having too large compartments shows as *excessive damping* (or *underestimation*) with regard to spatially-sensitive effects: namely spike propagation, and the contribution of point processes like synapses and probes to the neuron's state.
* On the other hand, having too small compartments may stress the numerical methods, over- or undershooting due to *numerical round-off* when thedifferences become tiny.

Generally, simulators don't like handling multiple timescales; keep things as rough as they still let your model work as designed.

..
    NEXT
    how do i different discretisations?? 
    --------------------------------------------

    let's port a model from a non-neuron first i guess.
    note that decisions are model specific...


How do I add calcium-modulated ion channels?
--------------------------------------------

With a custom `LEMS component <intro_lems.ipynb#Example:-Ca2⁺-dependent-base-conductance>`__.

How do I add Ornstein-Uhlenbeck noise?
--------------------------------------

With a custom `LEMS component <intro_lems.ipynb#Example:-Ornstein-Uhlenbeck-noise>`__.

..
    How do I add synapse timing dependent plasticity?
    -------------------------------------------------

    With a custom `LEMS component <example_stdp.ipynb>`__. NEXT


How do i add *other ion*-modulated ion channels or other mechanisms (modulated by possibly other quantities)?
-------------------------------------------------------------------------------------------------------------

With EDEN's LEMS extension `\<𝚅𝚊𝚛𝚒𝚊𝚋𝚕𝚎𝚁𝚎𝚏𝚎𝚛𝚎𝚗𝚌𝚎\> <extension_pointers.ipynb>`__.

How do I control mechanisms with a single "global" variable?
------------------------------------------------------------

With `<𝚅𝚊𝚛𝚒𝚊𝚋𝚕𝚎𝚁𝚎𝚏𝚎𝚛𝚎𝚗𝚌𝚎>s <extension_pointers.ipynb>`__ that point to the same single quantity (that can be on an abstract cell or wherever).

.. 
    LATER network element?

How do I explicitly set some variables on *every single timestep*?
------------------------------------------------------------------

By keeping a state variable ``timeSince`` and adding a tag like:

.. code-block:: xml
    
    <OnCondition test="t >= timeSince + (one timestep)">
        ... update things ...
        <StateAssignment variable="timeSince" value="t"/>
    </OnCondition>

See also the the `OU noise example <intro_lems.ipynb#Example:-Ornstein-Uhlenbeck-noise>`__ which updates the state variable ``i`` following a random walk.

EDEN and many other NeuroML-capable simulators also run ``<OnCondition test="(something unconditional)">`` on every timestep which saves one state variable, but it's not part of the specification.  Expect discrete updates to become part of the specification some day.


Where can I lean more about 


Writing LEMS equations
**********************

How do I add *normally distributed* random variables in my LEMS equations?
-----------------------------------------------------------------------------

A simple way is the `Box-Muller transform <https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform>`_.  Refer to the `OU noise example <intro_lems.ipynb#Example:-Ornstein-Uhlenbeck-noise>`__ for details.


How do I use the `sgn(x) <https://en.wikipedia.org/wiki/Sign_function>`_ function in expressions?
-------------------------------------------------------------------------------------------------------------

Use :math:`2*H(x) - 1` instead, and `request the feature <https://github.com/NeuroML/NeuroML2/issues>`_ if it would help a lot to have it.


What does ``log`` mean in equations?
------------------------------------

It stands for the *natural logarithm* (and is an alias for ``ln``), so that :math:`log(e) = 1`.  If you prefer the decimal logarithm, ``log10`` can be used instead.

~~~~

Troubleshooting
***************


Why is my spatial cell behaving exactly as if it was a single compartment?
--------------------------------------------------------------------------

Either:

* the cytoplasmic ``<resistivity>`` ``value`` is too low;
* the lengths of neurites in the ``Morphology`` are too small;
* the membrane's ``<specificCapacitance>`` is too low;
* the channels' conductance is too high;
* or (a cable in) the whole neuron was specified to be a single compartment with the `"unbranched section" <intro_spatial.ipynb#The-unbranched-section-directive>`__ directive, inadvertently.

My NeuroML model behaves very differently on EDEN compared to when it runs on e.g. NEURON via jNeuroML !
--------------------------------------------------------------------------------------------------------

Please `contact us <contact_us.rst>`__ to track down the issue and provide fixes and options.



When I use ``<`` in LEMS conditions, I get XML errors!
------------------------------------------------------

Unfortunately this character annoys XML. Use the `FORTRAN spelling <intro_lems.ipynb#The-components-of-\<Dynamics\>>`__ ``.lt.`` (less than) or ``.le.`` (less/equal) instead, or escape the character with ``&lt;`` (ugly).

~~~~

What is Eden capable of?
************************

Can EDEN run Pong?
------------------

Yes!  Check out the `example <example_pong.ipynb>`__. 

..
    Can EDEN run `D00M`?
    --------------------

    Perhaps in theory, but it won't be practical.  Instead you may run it as an OS process[capture doom?], and [feed video] to a neural network simulated in EDEN, to interact with the game. LATER

Does EDEN support Blockchain?
-----------------------------

In NeuroML, `anatomically detailed cells <intro_spatial.ipynb>`_ comprise spans of neurite.  Each such span can be considered as a chain of interlinked "blocks", also called "compartments" in our terminology.  In that sense, EDEN can simulate multiple blockchains within even a single neuron.

Does EDEN support Artificial Intelligence?
------------------------------------------

Artificial Intelligence (A.I.) may be implemented in various ways; one that's especially bio-inspired is *spiking neural networks*, whose activity can be simulated *in silico* using EDEN.
Since EDEN is primarily a scientific tool, it is up to the modeller to design effective A.I. architectures.

Can EDEN run Large Language Models?
-----------------------------------

.. 
    refer to ... LATER

As long they can be described as spiking neural networks, they can run.  The attention mechanism is not quite SNN-like, but it may be implementable in LEMS.  An example is in development.


What does Eden *not* do, and what are the workarounds?
******************************************************


Recording ``<DerivedVariable>``\ s
----------------------------------

Instead, add a shadow ``<StateVariable>`` for each ``<DerivedVariable>``; then for each shadow state variable, add a ``<OnCondition test="1 == 1">`` with a  ``<StateAssignment>`` of the shadow variable to the ``DerivedVariable``\ . (Recording ``<DerivedVariable>``\ s will be supported soon.)

``<Regime>``\ s
----------------

Instead, a ``<StateVariable>`` can be used to indicate the type of the regime, with the affected ``<DerivedVariable>``\ s becoming ``<ConditionalDerivedVariable>``\ s.

``<Display>`` elements in ``<Simulation>``
--------------------------------------------

Instead, use ``<OutputFile>`` and display with your preferred graphing system, as shown in the guide.

``<doubleSynapse>``\ s and ``<compositeInput>``\ s
--------------------------------------------------

Instead, merge the involved ``<ComponentType>``\ s into a single LEMS component that includes the dynamics of all parts, and exposes their summed influence.
