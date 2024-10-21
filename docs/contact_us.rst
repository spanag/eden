.. _contact_us:
    
##########
Contact us
##########

If you are interested in EDEN, we'll be delighted to hear from you:

.. table...

.. LATER forum
.. NEXT substitute chatroom

.. grid:: 1 1 2 2

    .. grid-item-card::
        :padding: 2
        :columns: 6

        **General discussion**
        ^^^
        Get in touch with our community on our chatroom and `mailing list <https://groups.google.com/g/eden-simulator-general>`_ - or just :ref:`email the developers <email_us>`.

    .. grid-item-card::
        :padding: 2
        :columns: 6

        **Support tickets**
        ^^^
        If you'd like to ask for a new feature, or report a problem, please file a support ticket on `GitLab <https://gitlab.com/c7859/neurocomputing-lab/Inferior_OliveEMC/eden/-/issues>`_.

.. https://stackoverflow.com/questions/13039131/applying-css-and-roles-for-text-blocks-instead-of-inline-spans-in-sphinx
.. rst-class:: sidenote
    
    To subscribe to the mailing list, send an email to eden-simulator-general+subscribe@googlegroups.com\ .


.. _email_us:

:hidden:`Email us`
******************

.. LATER invisible header here?

And more generally, feel free to send e-mail to the developers anytime:

    - Primary developer: Sotirios Panagiotou, s.panagiotou@erasmusmc.nl
    - Principal investigator: Dr.Christos Strydis, c.strydis@erasmusmc.nl

.. tip::
    
    For questions and suggestions regarding NeuroML, contact also the wider `NeuroML community <https://docs.neuroml.org/NeuroMLOrg/CommunicationChannels.html>`_.


Join us
*******

Do you want to add something to EDEN, like a feature, a modified engine or an improved user guide? Join us!

.. rst-class:: sidenote
    
    (See also: We have open :ref:`student thesis topics <thesis_topics>` on sim tech research!)

Guide contributions
===================

Do you have a nice usage example, or method of experimenting, with EDEN?
Or perhaps you'd like to see some part of this guide explained in more detail? |Contact us| for suggestions.

.. hint::
    
    Refer to :ref:`gallery` to see if something is missing.


Code contributions
==================

Would you like to add a new simulation feature, customize the existing codebase to your needs, or contrubute some useful API extensions or language bindings?  

See the Hacker's guide (TBA), |contact us| for support, and make a `pull request <https://gitlab.com/c7859/neurocomputing-lab/Inferior_OliveEMC/eden/-/merge_requests>`_.

.. _thesis_topics:

Research topics
===============

Neural simulation is a topic of active scientific research.  More than a product, EDEN is also designed as a platform to test new ideas on the technology of numerical simulation.  This ideas can range from new simulation algorithms, to adaptive performance tuning, to specialized computational hardware.

We do have several cutting-edge research questions on simulation technology; if a prospective student (or intern) has some relevant knowledge, and support from an academic institution or an intense personal interest, we can provide specific research topics to work on, and additional guidance.

.. See also the :ref:`hackers_guide` to understand the program's internals.
.. NEXT add topics for the general modelling workflow: model generation, results analysis...

Here are a few challenges, to explain our kind of work to prospective students:

- Get the sim code for a neural mass model from `GitHub <https://github.com/AmirrezaMov/tvb-algo-c/blob/master/tvb-algo.cpp>`_, and convert it so that the matrix is laid out in the `ELLPACK <https://faculty.cc.gatech.edu/~echow/ipcc/hpc-course/sparsemat.pdf>`_ format, for more than 1000 neurons and connection probabilities of 1%, 10%, 90%.  Compare simulation speed vs. the original code,  and and try to explain the results.
- Get a neural network model in the `ModelDB <https://modeldb.science/search?modeltype=Connectionist+Network%3B+Realistic+Network%3B+>`_ with more than, say, a thousand neurons.  Extract the set of synaptic connections between neurons.  Apply clustering algorithms to the neurons of the network, so that for N clusters the number of x * (max neurons per cluster) + y * (sum connections between clusters) is minimized.  Explore for varying x, y, N.

If you have pulled off one of these, or have another research project to propose, |contact us|; see also our `Msc thesis vacancies <https://neurocomputinglab.com/jobs/?job__type_spec=student-thesis>`_ (theme "BrainFrame"). 

.. NEXT add extensions? and backends... and benchmarks.
.. NEXT add refs to individual thesis topics?

