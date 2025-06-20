Pharmacodynamics
=================

While the mechanism of the pharmacokinetic is still poorly understood, the mechanism of actions of drugs at the molecular level is better understood (Bailey2005_). However, the link between the molecular level and the measured physiological variables is complex and thus, the pharmacodynamics models are usually empirical. The most common approach is to use the Hill function to describe the effect of the drug concentration on the physiological variables. The Hill function is a sigmoid function defined by the following equation:

.. math::

    E(t) = E_{max} \frac{C(t)^\gamma}{C(t)^\gamma + EC_{50}^\gamma}

where :math:`E(t)` is the effect of the drug at time :math:`t`, :math:`E_{max}` is the maximal effect, :math:`C(t)` is the drug concentration at time :math:`t`, :math:`\gamma` is the Hill coefficient, and :math:`EC_{50}` is the half-effect concentration (i.e., the concentration to obtain half the effect of the drugs). Figure below is an illustration of the sigmoid function for different values of the Hill coefficient and half-effect concentration.

.. figure:: ../images/sigmoid.png
   :width: 60%
   :align: center
   :alt: Sigmoid function

For propofol and remifentanil, before applying the Hill function, an effect-site compartment is added to the PK model to represent a delay between a rise of drug concentration in blood and the occurrence of the effect. This delay is dependent on the physiological variables, and thus, multiple effect-site compartments can be added to the model. As those compartments are virtual, the drug transfer is considered in only one direction, from blood to the effect site without affecting the blood compartment concentration. Thus, the addition of the effect site does not affect the PK model. The equation for one effect-site compartment is given by:

.. math::

    \dot{x}_{es}(t) = k_{e0} (x_1(t) - x_{es}(t))

where :math:`x_{es}(t)` is the drug concentration in the effect site, :math:`x_1(t)` is the drug concentration in blood, and :math:`k_{e0}` is the rate of drug transfer from blood to effect site. This leads to the full compartment model given in the next figure.

.. figure:: ../images/4_comportment_model.png
   :width: 70%
   :align: center
   :alt: Four-compartment model

   Four-compartment model for propofol and remifentanil.

In the simulator, we slightly abuse the notation and included the effect-site compartments in the PK model in order to keep all the dynamical system in the same state-space representation.  

BIS
------

If pharmacokinetics models usually assume no interaction between drugs, pharmacodynamics models should express the synergy or the antagonism between drugs. For the effect of propofol and remifentanil on the BIS, a 3D-Hill function is used to express the drug's synergy:

.. math::
    :label: eq:3DHill

    BIS(t) = BIS_{0} - E_{max} \frac{I(t)^\gamma}{1 + I(t)^\gamma}

with :math:`BIS_0` the initial BIS, :math:`E_{max}` the maximum effect of combined drugs, :math:`\gamma` the slope coefficient of the Hill curve and :math:`I(t)` the interaction term defined by:

.. math::

    I(t) = \frac{I_p(t) + I_r(t)}{1 - \beta \theta(t) + \beta \theta(t)^2}

where:

.. math::

    I_p(t) = \frac{x_{ep,BIS}(t)}{C_{50p,BIS}};\quad
    I_r(t) = \frac{x_{er,BIS}(t)}{C_{50r,BIS}};\quad
    \theta(t) = \frac{I_p(t)}{I_p(t)+I_r(t)}

In those equations, :math:`x_{ep,BIS}` and :math:`x_{er,BIS}` are the propofol and remifentanil concentrations of the BIS effect site, :math:`C_{50p,BIS}` and :math:`C_{50r,BIS}` are the propofol and remifentanil half-effect concentrations for BIS, and :math:`\beta` is the interaction term between the two drugs.

Few studies have been conducted on the pharmacodynamic part of the anesthesia process, and the models are less standardized. In this simulator, the values of the parameters of the 3D-Hill function are taken from the study of Bouillon2004_. The surface of the 3D-Hill function with the values from the mentioned study is shown in the figure below.

.. figure:: ../images/3Dhill.png
   :width: 80%
   :align: center
   :alt: 3D-Hill function


Tolerance of Laryngoscopy
-----------------------------

To output an indicator of analgesia in the simulator, we used the Tolerance of Laryngoscopy (TOL). The TOL is defined as the probability of reaction of the patient to the laryngoscopy. In Bouillon2004_, the authors proposed a hierarchical model to link drug effect site concentration to TOL. The model is given by:

.. math:: postopioid(t) = preopioid \times \left(1 - \frac{x_{er,BIS}(t)^{\gamma_r}}{x_{er,BIS}(t)^{\gamma_r} + (C_{r,50,TOL} \times preopioid)^{\gamma_r}}\right)
.. math:: TOL(t) = \frac{x_{ep,BIS}(t)^{\gamma_p}}{x_{ep,BIS}(t)^{\gamma_p} + (C_{p,50,TOL} \times postopioid(t))^{\gamma_p}}

where :math:`preopioid` is the tolerance of laryngoscopy without remifentanil, :math:`x_{er,BIS}(t)` and :math:`x_{ep,BIS}(t)` are the remifentanil and propofol concentration in the TOL effect site (same than the BIS effect site), :math:`C_{r,50,TOL}` and :math:`C_{p,50,TOL}` are the remifentanil and propofol half-effect concentrations for TOL.

Haemodynamic
--------------

For the effect of propofol and remifentanil on mean arterial pressure (MAP), the interaction of drugs has still to be studied. Thus, the effect of propofol, remifentanil and norepinephrine is considered to be independent and additive. The influence of propofol on MAP has been studied in Jeleazcov2015_, the influence of remifentanil in Standing2010_ and the one of norepinephrine in Beloeil2005_. For propofol, the authors of Jeleazcov2015_ find that the use of two different effect-site compartments better represents the effect of propofol on MAP. The model is given by:

.. math::
    \small
    MAP(t) =  MAP_0 - \underbrace{E_{max,r}\frac{x_{er,hemo}(t)^{\gamma_{r}}}{C_{50r,MAP}^{\gamma_{r}} + x_{er,hemo}(t)^{\gamma_{r}}}}_{\text{remifentanil effect}} 
    - \underbrace{E_{max,p}  \frac{I_p(t)}{1 + I_p(t)}}_{\text{propofol effect}} + \underbrace{E_{max,n}\frac{x_{n}(t)^{\gamma_{r}}}{C_{50n,MAP}^{\gamma_{n}} + x_{n}(t)^{\gamma_{r}}}}_{\text{norepinephrine effect}}

with:

.. math::

    I_p(t) = \left( \frac{x_{ep,hemo,1}(t)}{C_{50p,MAP,1}}\right)^{\gamma_{p1}} + \left(\frac{x_{ep,hemo,2}(t)}{C_{50p,MAP,2}}\right)^{\gamma_{p2}}

where :math:`MAP_0` is the MAP baseline, :math:`E_{max,r}`, :math:`E_{max,p}` and :math:`E_{max,n}` are the maximal effects of remifentanil, propofol and norepinephrine on MAP, :math:`x_{er,hemo}`, :math:`x_{ep,hemo,1}`, :math:`x_{ep,hemo,2}` and :math:`x_{n}` are the remifentanil and propofol, and norepinephrine concentrations in the hemodynamic effect site, or blood compartment for norepinephrine. :math:`C_{50r,MAP}`, :math:`C_{50p,MAP,1}`, :math:`C_{50p,MAP,2}`, :math:`C_{50n,MAP}`, :math:`\gamma_{r}`, :math:`\gamma_{p1}`, :math:`\gamma_{p2}`, and :math:`\gamma_{n}`, are the half-effect concentrations and Hill coefficients of remifentanil and propofol and norepinephrine.

For the effect on cardiac output (CO), studies are scarce. As for MAP, we considered additive drug effect, without any synergic effect. Because no sigmoid model was available in the litterature, we infer value to match experimental values from the following papers: Fairfield1991_ for propofol, Chanavaz2005_ for remifentanil and Monnet2011_ for  norepinephrine. We used the same effect sites than the one from MAP, and for propofol the mean concentration between the two efefct site compartment is used. Note that this is a crude simplification.

.. math::
    \small
    CO(t) =  CO_0 - \underbrace{E_{max,r}\frac{x_{er,hemo}(t)^{\gamma_{r}}}{C_{50r,CO}^{\gamma_{r}} + x_{er,hemo}(t)^{\gamma_{r}}}}_{\text{remifentanil effect}} 
    - \underbrace{E_{max,p}\frac{x_{p,hemo}(t)^{\gamma_{p}}}{C_{50p,CO}^{\gamma_{p}} + x_{p,hemo}(t)^{\gamma_{p}}}}_{\text{propofol effect}} + \underbrace{E_{max,n}\frac{x_{n}(t)^{\gamma_{r}}}{C_{50n,CO}^{\gamma_{n}} + x_{n}(t)^{\gamma_{r}}}}_{\text{norepinephrine effect}}

with:

.. math::

    x_{p,hemo}(t) = \frac{x_{ep,hemo,1}(t)+x_{ep,hemo,2}(t)}{2}

where :math:`CO_0` is the CO baseline, :math:`E_{max,r}`, :math:`E_{max,p}` and :math:`E_{max,n}` are the maximal effects of remifentanil, propofol and norepinephrine on CO, :math:`x_{er,hemo}`, :math:`x_{ep,hemo,1}`, :math:`x_{ep,hemo,2}` and :math:`x_{n}` are the remifentanil and propofol, and norepinephrine concentrations in the hemodynamic effect site, or blood compartment for norepinephrine. :math:`C_{50r,CO}`, :math:`C_{50p,CO}`, :math:`C_{50n,CO}`, :math:`\gamma_{r}`, :math:`\gamma_{p}`, and :math:`\gamma_{n}`, are the half-effect concentrations and Hill coefficients of remifentanil and propofol and norepinephrine.

The overall model of the anesthesia process is then given by connecting the PK model and the PD model. This can be formalized as a model with a linear dynamic and a non-linear output function in the following state-space representation:

.. math::
    :label: eq:standard_model

    \begin{cases}
        \dot{x}(t) = A x(t) + B u(t) \\
        y(t) = h(x(t))
    \end{cases}

where :math:`x(t)` is the system state, including the drug concentrations of propofol, remifentanil and norepinephrine in each compartment, :math:`u(t)` the drugs rates, and :math:`y(t)` is the output of the system, *i.e.*, the BIS, TOL, MAP, and CO.

Effect summary
-----------------

The following table summarizes the effect of the drugs on the physiological variables:

.. raw:: html

    <style>
      .blue-bg { background-color: #cce5ff; }  /* Light blue */
      .red-bg { background-color: #f8d7da; }   /* Light red */
      table.colored-table {
        border-collapse: collapse;
        width: 100%;
      }
      table.colored-table th,
      table.colored-table td {
        border: 1px solid #ddd;
        padding: 8px;
        text-align: center;
      }
      table.colored-table th {
        background-color: #f2f2f2;
      }
    </style>

    <table class="colored-table">
      <thead>
        <tr>
          <th></th>
          <th>Propofol</th>
          <th>Remifentanil</th>
          <th>Norepinephrine</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>BIS</th>
          <td class="blue-bg">-</td>
          <td class="blue-bg">-</td>
          <td>No effect</td>
        </tr>
        <tr>
          <th>TOL</th>
          <td class="red-bg">+</td>
          <td class="red-bg">+</td>
          <td>No effect</td>
        </tr>
        <tr>
          <th>MAP</th>
          <td class="blue-bg">-</td>
          <td class="blue-bg">-</td>
          <td class="red-bg">+</td>
        </tr>
        <tr>
          <th>CO</th>
          <td class="blue-bg">-</td>
          <td class="blue-bg">-</td>
          <td class="red-bg">+</td>
        </tr>
      </tbody>
    </table>

References
----------
.. [Bailey2005]     J. M. Bailey and W. M. Haddad, “Drug dosing control in clinical pharmacology,” *IEEE Control Systems Magazine*,
    vol. 25, no. 2, pp. 35–51, Apr. 2005, doi: https://doi.org/10.1109/MCS.2005.1411383.
.. [Bouillon2004] T. W. Bouillon et al., “Pharmacodynamic Interaction between Propofol and Remifentanil
    Regarding Hypnosis, Tolerance of Laryngoscopy, Bispectral Index, and Electroencephalographic Approximate
    Entropy,” Anesthesiology, vol. 100, no. 6, pp. 1353–1372, Jun. 2004, doi: https://doi.org/10.1097/00000542-200406000-00006.
.. [Ionescu2021] Ionescu, C. M., Neckebroek, M., Ghita, M., & Copot, D. (2021). An Open Source Patient
    Simulator for Design and Evaluation of Computer Based Multiple Drug Dosing Control for Anesthetic and
    Hemodynamic Variables. IEEE Access, 9, 8680–8694. https://doi.org/10.1109/ACCESS.2021.3049880
.. [Jeleazcov2015] C. Jeleazcov, M. Lavielle, J. Schüttler, and H. Ihmsen, “Pharmacodynamic response modelling
    of arterial blood pressure in adult volunteers during propofol anaesthesia,” BJA: British Journal of Anaesthesia,
    vol. 115, no. 2, pp. 213–226, Aug. 2015, doi: https://doi.org/10.1093/bja/aeu553.
.. [Standing2010] J. F. Standing, G. B. Hammer, W. J. Sam, and D. R. Drover, “Pharmacokinetic–pharmacodynamic
    modeling of the hypotensive effect of remifentanil in infants undergoing cranioplasty,” Pediatric Anesthesia,
    vol. 20, no. 1, pp. 7–18, 2010, doi: https://doi.org/10.1111/j.1460-9592.2009.03174.x.
..  [Beloeil2005]  H. Beloeil, J.-X. Mazoit, D. Benhamou, and J. Duranteau, “Norepinephrine kinetics and dynamics
    in septic shock and trauma patients,” BJA: British Journal of Anaesthesia, vol. 95, no. 6,
    pp. 782–788, Dec. 2005, doi: https://doi.org/10.Beloeil20051093/bja/aei259.
.. [Fairfield1991] J. E. Fairfield, A. Dritsas, and R. J. Beale, “Haemodynamic effects of propofol:
    induction with 2.5 mg/kg,” British Journal of Anaesthesia, vol. 67, no. 5,
    pp. 618–620, Nov. 1991, doi: https://doi.org/10.1093/bja/67.5.618.
.. [Chanavaz2005] C. Chanavaz et al., “Haemodynamic effects of remifentanil in children
    with and without intravenous atropine. An echocardiographic study,”
    BJA: British Journal of Anaesthesia, vol. 94, no. 1, pp. 74–79, Jan. 2005, doi: https://doi.org/10.1093/bja/aeh293.
..  [Monnet2011]  X. Monnet, J. Jabot, J. Maizel, C. Richard, and J.-L. Teboul, “Norepinephrine increases
    cardiac preload and reduces preload dependency assessed by passive leg raising in septic shock patients”
    Critical Care Medicine, vol. 39, no. 4, p. 689, Apr. 2011, doi: https://doi.org/10.1097/CCM.0b013e318206d2a3.
