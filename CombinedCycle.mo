package CombinedCycle
  package Substances
    model Gas
      import CombinedCycle.Types.*;
      input Pressure p(start = pstart);
      input Temperature T(start = Tstart);
      GasDensity d(start = dstart);
      SpecificEnthalpyGas h(start = hstart);
      SpecificEnergyGas u(start = u_start);
      constant Real R = Modelica.Constants.R;
      parameter Pressure pstart = 101325;
      parameter Temperature Tstart = 300;
      parameter GasDensity dstart = pstart * MM / Tstart / R;
      parameter SpecificEnthalpyGas hstart = cp * (Tstart - T0);
      parameter SpecificEnergyGas u_start = cv * (Tstart - T0);
      parameter Real MM = 0.029;
      parameter Temperature T0 = 298;
      parameter SpecificHeatCapacity cp = 1000;
      final parameter SpecificHeatCapacity cv = cp - R / MM;
    equation
      p / d = R / MM * T;
      h / hstart = cp * (T - T0) / hstart;
      u / u_start = cv * (T - T0) / u_start;
      annotation(Documentation(info = "<HTML>
<p>Gas model with Ideal Gas state equation and constant specific heat capacity</p>
</HTML>"));
    end Gas;

    model WaterLiquid
      import CombinedCycle.Types.*;
      //Reference parameters for adimensionalization
      constant Pressure pcrit = 22.064e6;
      constant Temperature Tcrit = 647.096;
      constant SpecificHeatCapacity R = 461.522;
      constant SpecificEnthalpyLiquid href = R * Tcrit;
      constant LiquidDensity dref = pcrit / (R * Tcrit);
      //Start parameters
      parameter Pressure pstart = 1.5e5;
      parameter Temperature Tstart = 350;
      parameter Real pr_start = pstart / pcrit;
      parameter Real Tr_start = Tstart / Tcrit;
      parameter Real pr_sat_start = ((0.803088 * exp(Tr_start) + (-3.96047)) * exp(Tr_start) + 6.51728) * exp(Tr_start) + (-3.57666);
      parameter Real h_lsr_start = 5.55738 + log(pr_sat_start) * (1.35186 + log(pr_sat_start) * (0.142999 + log(pr_sat_start) * 0.00646981));
      parameter Real d_lsr_start = 7.88128 + log(pr_sat_start) * ((-2.39458) + log(pr_sat_start) * ((-0.477231) + log(pr_sat_start) * ((-0.0500329) + log(pr_sat_start) * (-0.00211171))));
      parameter Pressure p_sat_start = pr_start * pcrit;
      parameter SpecificEnthalpyLiquid h_ls_start = h_lsr_start * href;
      parameter LiquidDensity d_start = d_lsr_start * dref;
      parameter SpecificEnthalpyLiquid h_start = h_ls_start + (pstart - p_sat_start) / d_start;
      parameter SpecificEnergyLiquid u_start = h_start - pstart / d_start;
      //Variables
      input Pressure p(start = pstart);
      input Temperature T(start = Tstart);
      Real Tr(start = Tr_start);
      Real pr_sat(start = pr_sat_start);
      Real h_lsr(start = h_lsr_start);
      Real d_lsr(start = d_lsr_start);
      Pressure p_sat(start = p_sat_start);
      SpecificEnthalpyLiquid h_ls(start = h_ls_start);
      LiquidDensity d(start = d_lsr_start);
      SpecificEnthalpyLiquid h(start = h_start);
      SpecificEnergyLiquid u(start = u_start);
    equation
      Tr = T / Tcrit;
      pr_sat = ((0.803088 * exp(Tr) + (-3.96047)) * exp(Tr) + 6.51728) * exp(Tr) + (-3.57666);
      //Tr = 0.989276 + log(pr_sat) * (0.123086 + log(pr_sat) * (0.0108249 + log(pr_sat) * 0.000405073));
      h_lsr = 5.51886 + log(pr_sat) * (1.31355 + log(pr_sat) * (0.131269 + log(pr_sat) * 0.0053645));
      d_lsr = 7.93122 + log(pr_sat) * ((-2.32899) + log(pr_sat) * ((-0.446718) + log(pr_sat) * ((-0.0440744) + log(pr_sat) * (-0.00169849))));
      p_sat = pr_sat * pcrit;
      h_ls = h_lsr * href;
      d = d_lsr * dref;
      h = h_ls + (p - p_sat) / d;
      u = h_ls - p_sat / d;
      assert(p >= p_sat, "This model is valid only if p>=p saturation.\nSaturation pressure is: " + String(p_sat) + "[Pa], while the pressure is: " + String(p) + "[Pa].");
      annotation(Documentation(info = "<HTML>
<p> Liquid water model, for a given temperature the saturated liquid state properties are calculated and then corrected in regard to pressure treating the fluid if it was incompressible.</p>
<p> The functions are obtained by linear regression of the IF97 model.</p>
</HTML>"));
    end WaterLiquid;

    model WaterSaturation
      import CombinedCycle.Types.*;
      //Reference parameters for adimensionalization
      constant Pressure pcrit = 22.064e6;
      constant Pressure Tcrit = 647.096;
      constant SpecificHeatCapacity R = 461.522;
      constant SpecificEnthalpyGas href = R * Tcrit;
      constant GasDensity dref = pcrit / (R * Tcrit);
      //Saturation properties and relative reduced properties
      input Pressure p(start = pstart);
      Real pr(start = pr_start) "Reduced Temperature";
      Temperature T_s(start = T_s_start) "Saturation temperature";
      Real Tr(start = Tr_start) "Reduced Temperature";
      LiquidDensity d_ls(start = d_ls_start) "Density of liquid at saturation";
      Real d_lsr(start = d_lsr_start) "Reduced density of liquid at saturation";
      GasDensity d_vs(start = d_vs_start) "Density of vapour at saturation";
      Real d_vsr(start = d_vsr_start) "Reduced density of vapour at saturation";
      SpecificEnthalpyLiquid h_ls(start = h_ls_start);
      Real h_lsr(start = h_lsr_start);
      SpecificEnthalpyGas h_vs(start = h_vs_start);
      Real h_vsr(start = h_vsr_start);
      VolumetricEnergyLiquid de_ls(start = de_ls_start, nominal = 1e9) "Density-specific energy product for liquid at saturation condition";
      Real de_lsr(start = de_lsr_start, nominal = 1e9 / pcrit) "Reduced density-specific energy product for liquid at saturation condition";
      VolumetricEnergyGas de_vs(start = de_vs_start, nominal = 1e7) "Density-specific energy product for vapour at saturation condition";
      Real de_vsr(start = de_vsr_start, nominal = 1e7 / pcrit) "Reduced density-specific energy product for vapour at saturation condition";
      DerSaturationTemperatureByPressure d_T_s_dp(start = d_T_s_dp_start, nominal = 1e-5);
      Real d_T_s_dp_r(start = d_T_s_dp_r_start, nominal = 1e-5 * pcrit / Tcrit);
      DerDensityByPressureLiquid d_d_ls_dp(start = d_d_ls_dp_start, nominal = 1e-3);
      Real d_d_ls_dp_r(start = d_d_ls_dp_r_start, nominal = 1e-3 * href);
      DerDensityByPressureGas d_d_vs_dp(start = d_d_vs_dp_start, nominal = 1e-6);
      Real d_d_vs_dp_r(start = d_d_vs_dp_r_start, nominal = 1e-6 * href);
      DerVolumetricEnergyByPressureLiquid d_de_ls_dp(start = d_de_ls_dp_start, nominal = 1e3);
      DerVolumetricEnergyByPressureGas d_de_vs_dp(start = d_de_vs_dp_start, nominal = 10);
      SpecificEnergyLiquid u_ls(start = u_ls_start, nominal = 1e6);
      SpecificEnergyGas u_vs(start = u_vs_start, nominal = 1e6);
      //Start parameters
      parameter Pressure pstart = 1.5e5;
      parameter Real pr_start = pstart / pcrit;
      parameter Temperature T_s_start = Tr_start * Tcrit;
      parameter LiquidDensity d_ls_start = d_lsr_start * dref;
      parameter GasDensity d_vs_start = d_vsr_start * dref;
      parameter SpecificEnthalpyLiquid h_ls_start = h_lsr_start * href;
      parameter SpecificEnthalpyGas h_vs_start = h_vsr_start * href;
      parameter VolumetricEnergyLiquid de_ls_start = de_lsr_start * pcrit;
      parameter VolumetricEnergyGas de_vs_start = de_vsr_start * pcrit;
      parameter DerSaturationTemperatureByPressure d_T_s_dp_start = d_T_s_dp_r_start / pcrit * Tcrit;
      parameter DerDensityByPressureLiquid d_d_ls_dp_start = d_d_ls_dp_r_start / href;
      parameter DerDensityByPressureGas d_d_vs_dp_start = d_d_vs_dp_r_start / href;
      parameter Real Tr_start = 0.991429 + log(pr_start) * (0.125239 + log(pr_start) * (0.0114884 + log(pr_start) * 0.000468091));
      parameter Real d_lsr_start = 7.88128 + log(pr_start) * ((-2.39458) + log(pr_start) * ((-0.477231) + log(pr_start) * ((-0.0500329) + log(pr_start) * (-0.00211171))));
      parameter Real d_vsr_start = 0.00230966 + pr_start * (1.46895 + pr_start * 0.0760944);
      parameter Real h_lsr_start = 5.55738 + log(pr_start) * (1.35186 + log(pr_start) * (0.142999 + log(pr_start) * 0.00646981));
      parameter Real h_vsr_start = 8.92734 + log(pr_start) * ((-0.559025) + log(pr_start) * ((-0.215448) + log(pr_start) * ((-0.0277462) + log(pr_start) * (-0.00126407))));
      parameter Real de_lsr_start = 50.4827 + log(pr_start) * (6.71714 + log(pr_start) * (0.0550618 + log(pr_start) * (-0.0139366)));
      parameter Real de_vsr_start = 0.00786636 + pr_start * (13.3487 + pr_start * ((-6.48991) + pr_start * 26.2566));
      parameter DerSaturationTemperatureByPressure d_T_s_dp_r_start = (0.125239 + log(pr_start) * (0.0229768 + log(pr_start) * 0.001404273)) / pr_start;
      parameter Real d_d_ls_dp_r_start = ((-2.39458) + log(pr_start) * ((-0.954462) + log(pr_start) * ((-0.1500987) + log(pr_start) * (-0.00844684)))) / pr_start;
      parameter Real d_d_vs_dp_r_start = 1.46895 + pr_start * 0.1521888;
      parameter DerVolumetricEnergyByPressureLiquid d_de_ls_dp_start = (6.71714 + log(pr_start) * (0.1101236 + log(pr_start) * (-0.0418098))) / pr_start;
      parameter DerVolumetricEnergyByPressureGas d_de_vs_dp_start = 13.3487 + pr_start * ((-12.97982) + pr_start * 78.7698);
      parameter SpecificEnergyLiquid u_ls_start = h_ls_start - pstart / d_ls_start;
      parameter SpecificEnergyGas u_vs_start = h_vs_start - pstart / d_vs_start;
    equation
      //Reduced Properties
      pr = p / pcrit;
      Tr = T_s / Tcrit;
      d_lsr = d_ls / dref;
      d_vsr = d_vs / dref;
      h_lsr = h_ls / href;
      h_vsr = h_vs / href;
      de_lsr = de_ls / pcrit;
      de_vsr = de_vs / pcrit;
      d_T_s_dp_r = d_T_s_dp * pcrit / Tcrit;
      d_d_ls_dp_r = d_d_ls_dp * href;
      d_d_vs_dp_r = d_d_vs_dp * href;
      //Correlations
      Tr = 0.989276 + log(pr) * (0.123086 + log(pr) * (0.0108249 + log(pr) * 0.000405073));
      d_lsr = 7.93122 + log(pr) * ((-2.32899) + log(pr) * ((-0.446718) + log(pr) * ((-0.0440744) + log(pr) * (-0.00169849))));
      d_vsr = 0.00227499 + pr * (1.46971 + pr * 0.0726054);
      h_lsr = 5.51886 + log(pr) * (1.31355 + log(pr) * (0.131269 + log(pr) * 0.0053645));
      h_vsr = 8.96044 + log(pr) * ((-0.515575) + log(pr) * ((-0.19525) + log(pr) * ((-0.0238049) + log(pr) * (-0.000990982))));
      de_lsr = 50.7317 + log(pr) * (6.96248 + log(pr) * (0.129285 + log(pr) * (-0.00704616)));
      de_vsr = 0.0076068 + pr * (13.3595 + pr * ((-6.60853) + pr * 26.636));
      d_T_s_dp_r = (0.123086 + log(pr) * (0.0216498 + log(pr) * 0.00121522)) / pr;
      d_d_ls_dp_r = ((-2.32899) + log(pr) * ((-0.893437) + log(pr) * ((-0.132223) + log(pr) * (-0.00679395)))) / pr;
      d_d_vs_dp_r = 1.46971 + pr * 0.145211;
      d_de_ls_dp = (6.96248 + log(pr) * (0.258569 + log(pr) * (-0.0211385))) / pr;
      d_de_vs_dp = 13.3595 + pr * ((-13.2171) + pr * 79.9081);
      u_ls = h_ls - p / d_ls;
      u_vs = h_vs - p / d_vs;
      annotation(Documentation(info = "<HTML>
<p> Liquid and vapour water at saturation properties as a function of pressure.
<p> The functions are obtained by linear regression of the IF97 model.
</HTML>"));
    end WaterSaturation;

    model WaterVapour
      import CombinedCycle.Types.*;
      //Reference parameters for adimensionalization
      constant Pressure pcrit = 22.064e6;
      constant Pressure Tcrit = 647.096;
      constant SpecificHeatCapacity R = 461.522;
      constant SpecificEnthalpyGas href = R * Tcrit;
      constant GasDensity dref = pcrit / (R * Tcrit);
      //Start parameters
      parameter Pressure pstart = 1.5e5;
      parameter Temperature Tstart = 400;
      parameter GasDensity d_start = dr_start * dref;
      parameter SpecificEnthalpyGas h_start = hr_start * href;
      parameter SpecificEnergyGas u_start = h_start - pstart / d_start;
      parameter Real pr_start = pstart / pcrit;
      parameter Real Tr_start = Tstart / Tcrit;
      parameter Real hr_start = hr_sat_start + Tr_start + (c1 * Tr_start + (a1 + 1) * c2 / Tr_start ^ (a1 - 1) + (a2 + 1) * c3 / Tr_start ^ (a2 - 1)) * dr_start + (c5 * Tr_start + (a3 + 2) * c6 / 2 / Tr_start ^ (a3 - 1) + (a4 + 2) * c7 / 2 / Tr_start ^ (a4 - 1)) * dr_start ^ 2 - (Tr_sat_start + (c1 * Tr_sat_start + (a1 + 1) * c2 / Tr_sat_start ^ (a1 - 1) + (a2 + 1) * c3 / Tr_sat_start ^ (a2 - 1)) * dr_start + (c5 * Tr_sat_start + (a3 + 2) * c6 / 2 / Tr_sat_start ^ (a3 - 1) + (a4 + 2) * c7 / 2 / Tr_sat_start ^ (a4 - 1)) * dr_start ^ 2) + (-2.8866624E-02 * (Tr_start ^ 4 - Tr_sat_start ^ 4)) + 2.3034065E-01 * (Tr_start ^ 3 - Tr_sat_start ^ 3) - 7.2944159E-02 * (Tr_start ^ 2 - Tr_sat_start ^ 2) + 2.9706549 * (Tr_start - Tr_sat_start);
      parameter Real dr_start = ((Tr_start ^ 2 + 4 * pr_start * (Tr_start * (c1 + c2 / Tr_start ^ a1 + c3 / Tr_start ^ a2 + c4 / Tr_start))) ^ 0.5 - Tr_start) / 2 / (Tr_start * (c1 + c2 / Tr_start ^ a1 + c3 / Tr_start ^ a2 + c4 / Tr_start));
      parameter Real hr_sat_start = 8.92734 + log(pr_sat_start) * ((-0.559025) + log(pr_sat_start) * ((-0.215448) + log(pr_sat_start) * ((-0.0277462) + log(pr_sat_start) * (-0.00126407))));
      parameter Real Tr_sat_start = 0.991429 + log(pr_sat_start) * (0.125239 + log(pr_sat_start) * (0.0114884 + log(pr_sat_start) * 0.000468091));
      parameter Real pr_sat_start = (-9.652156) + (93.1326 + dr_start / 0.0760944) ^ 0.5;
      //Constants of the correlations
      parameter Real a1 = 0.562;
      parameter Real a2 = 4.20;
      parameter Real a3 = 4.97;
      parameter Real a4 = 9.22;
      parameter Real c1 = 0.10052348;
      parameter Real c2 = 0.21589098;
      parameter Real c3 = -0.16400788;
      parameter Real c4 = -0.51121127;
      parameter Real c5 = 6.82878135e-2;
      parameter Real c6 = 0.100147264;
      parameter Real c7 = -2.44119951e-2;
      parameter Real c8 = -9.75556323e-2;
      //Variables
      input Pressure p(start = pstart);
      input Temperature T(start = Tstart);
      GasDensity d(start = d_start);
      SpecificEnthalpyGas h(start = h_start);
      SpecificEnergyGas u(start = u_start);
      Real pr(start = pr_start);
      Real Tr(start = Tr_start);
      Real hr(start = hr_start);
      Real dr(start = dr_start);
      Real hr_sat(start = hr_sat_start);
      Real Tr_sat(start = Tr_sat_start);
      Real pr_sat(start = pr_sat_start);
      //Real cp;
    equation
      //Adimensionalization
      pr = p / pcrit;
      Tr = T / Tcrit;
      hr = h / href;
      dr = d / dref;
      //Correlations
      pr = Tr * dr * (1 + (c1 + c2 / Tr ^ a1 + c3 / Tr ^ a2 + c4 / Tr) * dr + (c5 + c6 / Tr ^ a3 + c7 / Tr ^ a4 + c8 / Tr) * dr ^ 2);
      dr = 0.00227499 + pr_sat * (1.46971 + pr_sat * 0.0726054);
      Tr_sat = 0.989276 + log(pr_sat) * (0.123086 + log(pr_sat) * (0.0108249 + log(pr_sat) * 0.000405073));
      hr_sat = 8.96044 + log(pr_sat) * ((-0.515575) + log(pr_sat) * ((-0.19525) + log(pr_sat) * ((-0.0238049) + log(pr_sat) * (-0.000990982))));
      hr = hr_sat + Tr + (c1 * Tr + (a1 + 1) * c2 / Tr ^ (a1 - 1) + (a2 + 1) * c3 / Tr ^ (a2 - 1)) * dr + (c5 * Tr + (a3 + 2) * c6 / 2 / Tr ^ (a3 - 1) + (a4 + 2) * c7 / 2 / Tr ^ (a4 - 1)) * dr ^ 2 - (Tr_sat + (c1 * Tr_sat + (a1 + 1) * c2 / Tr_sat ^ (a1 - 1) + (a2 + 1) * c3 / Tr_sat ^ (a2 - 1)) * dr + (c5 * Tr_sat + (a3 + 2) * c6 / 2 / Tr_sat ^ (a3 - 1) + (a4 + 2) * c7 / 2 / Tr_sat ^ (a4 - 1)) * dr ^ 2) + (-2.8866624E-02 * (Tr ^ 4 - Tr_sat ^ 4)) + 2.3034065E-01 * (Tr ^ 3 - Tr_sat ^ 3) - 7.2944159E-02 * (Tr ^ 2 - Tr_sat ^ 2) + 2.9706549 * (Tr - Tr_sat);
      //cp = R * (1 + (c1 - (a1 ^ 2 - 1) * c2 / Tr ^ a1 - (a2 ^ 2 - 1) * c3 / Tr ^ a2) * dr + (c5 - (a3 + 2) * (a3 - 1) * c6 / 2 / Tr ^ a3 - (a4 + 2) * (a4 - 1) * c7 / 2 / Tr ^ a4) * dr ^ 2 + (-5.3290329e-2 * Tr ^ 3) + 0.3189217 * Tr ^ 2 - 6.7330668e-2 * Tr + 1.3710226);
      u = h - p / d;
      //assert(pr <= pr_sat, "This model is valid only if p<=p saturation.\n Reduced saturation pressure is: " + String(pr_sat) + ", while the reduced pressure is: " + String(pr) + ".");
      annotation(Documentation(info = "<HTML>
<p> Vapour water model, first the saturation properties are calculated using functions obtained by linear regression of the IF97 model. Then the actual properties are obtained with a cubic equation of state based on a trucnated virail expansion.
<p> Source: JLM Fernandes. Fast evaluation of thermodynamic properties of superheated steam: a cubic equation of state. Applied Thermal Engineering, 16(1):71-79, 1996.
</HTML>"));
    end WaterVapour;
  end Substances;

  package Connectors
    connector FlangeA "A-type flange connector for gas flows"
      Types.Pressure p(start = 1e5) "Pressure";
      flow Types.MassFlowRate w "Mass flowrate";
      input Types.SpecificEnthalpy h "Specific enthalpy of fluid";
      annotation(Icon(graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {159, 159, 223}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid)}), Documentation(info = "<HTML>
<p> Must always be connected to a single type-B connector <tt>FlangeB</tt>.
</HTML>", revisions = "<html>
<ul>
<li><i>20 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
    end FlangeA;

    connector FlangeB "B-type flange connector for gas flows"
      Types.Pressure p(start = 1e5) "Pressure";
      flow Types.MassFlowRate w "Mass flowrate";
      output Types.SpecificEnthalpy h "Specific enthalpy of fluid";
      annotation(Icon(graphics = {Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {159, 159, 223}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-42, 44}, {44, -40}}, lineColor = {159, 159, 223}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Documentation(info = "<HTML>
<p> Must always be connected to a single type-A connector <tt>FlangeA</tt>.
</HTML>", revisions = "<html>
<ul>
<li><i>20 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
    end FlangeB;

    connector HeatPort_a "Thermal port for 1-dim. heat transfer (filled rectangular icon)"
      Types.Temperature T "Port temperature";
      flow Types.HeatFlowRate Q_flow "Heat flow rate (positive if flowing from outside into the component)";
      annotation(defaultComponentName = "port_a", Documentation(info = "<HTML>
<p>This connector is used for 1-dimensional heat flow between components.
The variables in the connector are:</p>
<pre>   
   T       Temperature in [Kelvin].
   Q_flow  Heat flow rate in [Watt].
</pre>
<p>According to the Modelica sign convention, a <b>positive</b> heat flow
rate <b>Q_flow</b> is considered to flow <b>into</b> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <b>HeatPort_a</b> and
<b>HeatPort_b</b> are identical with the only exception of the different
<b>icon layout</b>.</p></HTML>
"), Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid)}), Diagram(graphics = {Rectangle(extent = {{-50, 50}, {50, -50}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-120, 120}, {100, 60}}, lineColor = {191, 0, 0}, textString = "%name")}));
    end HeatPort_a;

    connector HeatPort_b "Thermal port for 1-dim. heat transfer (unfilled rectangular icon)"
      Types.Temperature T "Port temperature";
      flow Types.HeatFlowRate Q_flow "Heat flow rate (positive if flowing from outside into the component)";
      annotation(defaultComponentName = "port_b", Documentation(info = "<HTML>
<p>This connector is used for 1-dimensional heat flow between components.
The variables in the connector are:</p>
<pre>
   T       Temperature in [Kelvin].
   Q_flow  Heat flow rate in [Watt].
</pre>
<p>According to the Modelica sign convention, a <b>positive</b> heat flow
rate <b>Q_flow</b> is considered to flow <b>into</b> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <b>HeatPort_a</b> and
<b>HeatPort_b</b> are identical with the only exception of the different
<b>icon layout</b>.</p></HTML>
"), Diagram(graphics = {Rectangle(extent = {{-50, 50}, {50, -50}}, lineColor = {191, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-100, 120}, {120, 60}}, lineColor = {191, 0, 0}, textString = "%name")}), Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {191, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}));
    end HeatPort_b;

    connector DHT "Distributed Heat Terminal"
      Types.Temperature T[N] "Temperature at the nodes";
      flow Types.HeatFlux phi[N] "Heat flux at the nodes";
      parameter Integer N(min = 1) = 2 "Number of nodes";
      annotation(Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {255, 127, 0}, fillColor = {255, 127, 0}, fillPattern = FillPattern.Solid)}));
    end DHT;

    connector MultiHeatPort
      parameter Integer N(min = 1) = 3;
      CombinedCycle.Connectors.HeatPort_a hp[N];
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(lineColor = {255, 238, 44}, fillColor = {255, 238, 44}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}), Documentation(info = "<HTML>
<p> Heat port for multiple temperature and heat exchanges nodes.
</HTML>"));
    end MultiHeatPort;

    connector HFP "Heat flux port"
      Types.Temperature T "Temperature at the nodes";
      flow Types.HeatFlux phi "Heat flux at the nodes";
      annotation(Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {255, 127, 0}, fillColor = {255, 127, 0}, fillPattern = FillPattern.Solid)}), Documentation(info = "<HTML>
<p>This connector is used for 1-dimensional heat flux between components.
The variables in the connector are:</p>
<pre>   
   T       Temperature in [Kelvin].
   phi  Heat flux in [Watt / m^2].
</pre>
<p>According to the Modelica sign convention, a <b>positive</b> heat flow
rate <b>Q_flow</b> is considered to flow <b>into</b> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<b>icon layout</b>.</p></HTML>
"));
    end HFP;
  end Connectors;

  package Components
    package Gas
      model SourceW "Flowrate source for gas flows"
        extends CombinedCycle.Icons.Gas.SourceW;
        parameter Types.Pressure pnom = 101325 "Nominal pressure";
        parameter Types.Temperature Tnom = 300 "Nominal temperature";
        parameter Types.MassFlowRate wnom = 0 "Nominal mass flowrate";
        parameter Real G = 1 "HydraulicConductance";
        Types.MassFlowRate w;
        CombinedCycle.Substances.Gas gas(pstart = pnom, Tstart = Tnom);
        CombinedCycle.Connectors.FlangeB flange(w(start = wnom)) annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
      equation
        flange.w = (-w) + (flange.p - pnom) * G;
        w = wnom "Flow rate set by parameter";
        gas.p = flange.p;
        gas.h = flange.h;
        gas.T = Tnom "Temperature set by parameter";
        annotation(Icon(graphics), Documentation(info = "<html>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the nominal gas composition is given by <tt>Xnom</tt>,whose default value is <tt>Medium.reference_X</tt> .
<p>If <tt>G</tt> is set to zero, the flowrate source is ideal; otherwise, the outgoing flowrate decreases proportionally to the outlet pressure.</p>
<p>If the <tt>in_w0</tt> connector is wired, then the source massflowrate is given by the corresponding signal, otherwise it is fixed to <tt>w0</tt>.</p>
<p>If the <tt>in_T</tt> connector is wired, then the source temperature is given by the corresponding signal, otherwise it is fixed to <tt>T</tt>.</p>
<p>If the <tt>in_X</tt> connector is wired, then the source massfraction is given by the corresponding signal, otherwise it is fixed to <tt>Xnom</tt>.</p>
</html>", revisions = "<html>
<ul>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>w0fix</tt> and <tt>Tfix</tt> and <tt>Xfix</tt>; the connection of external signals is now detected automatically.</li> <br> Adapted to Modelica.Media
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
      end SourceW;

      model SourceP "Pressure source for gas flows"
        extends CombinedCycle.Icons.Gas.SourceP;
        parameter Types.Pressure pnom = 101325 "Nominal pressure";
        parameter Real Res = 1 "Hydraulic resistance";
        parameter Types.Temperature Tnom = 300 "Nominal temperature";
        CombinedCycle.Substances.Gas gas(p = pnom, T = Tnom);
        CombinedCycle.Connectors.FlangeB flange annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
      equation
        flange.p = gas.p + flange.w * Res;
        flange.h = gas.h;
        annotation(Icon(graphics), Diagram(graphics), Documentation(info = "<html>
  <p>Mass flow source of the Substance Gas with fixed Temperature and Pressure</p>
  </html>"));
      end SourceP;

      model SinkP "Pressure sink for gas flows"
        extends CombinedCycle.Icons.Gas.SourceP;
        parameter Types.Pressure pnom = 101325 "Nominal pressure";
        parameter Real Res = 1 "Hydraulic Resistance";
        CombinedCycle.Connectors.FlangeA flange annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
      equation
        flange.p = pnom + flange.w * Res;
        annotation(Icon(graphics), DymolaStoredErrors, Diagram(graphics), Documentation(info = "<html>
<p>Mass flow sink of the Substance Gas with fixed Temperature and Pressure</p>
<p>If <tt>R</tt> is set to zero, the pressure sink is ideal; otherwise, the inlet pressure increases proportionally to the outgoing flowrate.</p>
</html>"));
      end SinkP;

      model Flow1DGasPConst "Model for gas flows - no mass storage"
        extends CombinedCycle.Icons.Gas.Tube;
        import Modelica.SIunits.*;
        parameter Pressure pstart = 101325 "Pressure start value of outgoing gas";
        parameter Temperature Tstart = 300 "Temperature start value of outgoing gas";
        parameter Temperature Tmstart = 300 "Metal wall start temperature";
        parameter Volume V = 1 "Inner volume";
        Modelica.SIunits.Mass M(start = 30 * V, nominal = 30 * V) "Gas total mass";
        CombinedCycle.Substances.Gas gas_in;
        CombinedCycle.Substances.Gas gas_out(pstart = pstart, Tstart = Tstart);
        Real der_E(nominal = 1e6);
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Connectors.HeatPort_a wall(T(start = Tstart)) annotation(Placement(transformation(extent = {{-20, 42}, {20, 82}}, rotation = 0), iconTransformation(extent = {{-40, 40}, {40, 60}})));
      equation
        //Definitions
        M = gas_out.d * V;
        der_E = der(gas_out.T) * (gas_out.cv * M);
        0 = inlet.w + outlet.w "Mass balance";
        der_E = inlet.w * inlet.h + outlet.w * outlet.h + wall.Q_flow "Energy balance";
        inlet.p = outlet.p "Momentum balance";
        // Boundary conditions
        inlet.p = gas_in.p;
        outlet.p = gas_out.p;
        inlet.h = gas_in.h;
        outlet.h = gas_out.h;
        wall.T = (gas_out.T + gas_in.T) / 2;
        annotation(Placement(transformation(extent = {{-60, 40}, {60, 60}}, rotation = 0)), Icon(graphics = {Rectangle(extent = {{-50, 4}, {32, 0}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, 12}, {32, -8}, {52, 2}, {32, 12}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Documentation(info = "<html>
<p>This model describes a constant volume buffer with metal walls. The metal wall temperature and the heat transfer coefficient between the wall and the fluid are uniform. The wall is thermally insulated from the outside.</p>
<p>If the inlet or the outlet are connected to a bank of tubes, the model can actually represent a collector or a distributor.</p>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the start composition is given by <tt>Xstart</tt>, whose default value is <tt>Medium.reference_X</tt> .
</html>", revisions = "<html>
<ul>
<li><i>30 May 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Initialisation support added.</li>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>
"), Diagram(graphics));
      end Flow1DGasPConst;

      model BraytonJoule_exhaustGases "Flowrate source for modeling Brayton Joule exhaust gases, depending on the usage percentage"
        extends CombinedCycle.Icons.Gas.SourceW;
        parameter Types.Pressure pstart = 101325 "Starting pressure";
        parameter Types.Temperature T_min "Temperature at zero load";
        parameter Types.Temperature T_max "Max temperature";
        parameter Types.MassFlowRate w_min "Minimum flow rate";
        parameter Types.MassFlowRate w_max "Mass Flow rate at full load";
        parameter Real load_min_w = 0.5 "Load factor at minimum flow rate";
        parameter Types.MassFlowRate wstart = 0 "Starting mass flowrate";
        parameter Types.Temperature Tstart = 300 "Starting temperature";
        parameter Real sf = 1e2 "Smoothing factor (higher->steeper)";
        constant Real pi = Modelica.Constants.pi;
        Types.MassFlowRate w;
        CombinedCycle.Substances.Gas gas(pstart = pstart, Tstart = Tstart);
        CombinedCycle.Connectors.FlangeB flange(w(start = wstart, nominal = w_min)) annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput load "Percentage of usage" annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        Real lambda(start = 1) "One above load_min_w, zero below load_min_w, with smoothing";
      equation
        lambda = (atan(sf * (load - load_min_w)) + pi / 2) / pi;
        flange.w = -w;
        w = w_min * (1 - lambda) + (w_min + (w_max - w_min) * (load - load_min_w) * 2) * lambda;
        gas.T = (T_min + (T_max - T_min) * (load * 2)) * (1 - lambda) + T_max * lambda;
        flange.p = gas.p;
        flange.h = gas.h;
        annotation(Icon(graphics), Documentation(info = "<html>
<p>Flowrate source for modeling Brayton Joule exhaust gases. As the load increases the mass flow decreses and the Gas temperature increases following a smoothed out step function.
</html>"), Diagram(graphics));
      end BraytonJoule_exhaustGases;

      model BraytonJoule_exhaustGases_disturbances "Flowrate source for modeling Brayton Joule exhaust gases, depending on the usage percentage, with disturbances on temperature and flow rate"
        extends CombinedCycle.Icons.Gas.SourceW;
        parameter Types.Pressure pstart = 101325 "Starting pressure";
        parameter Types.Temperature T_min "Temperature at zero load";
        parameter Types.Temperature T_max "Max temperature";
        parameter Types.MassFlowRate w_min "Minimum flow rate";
        parameter Types.MassFlowRate w_max "Mass Flow rate at full load";
        parameter Real load_min_w = 0.5 "Load factor at minimum flow rate";
        parameter Types.MassFlowRate wstart = 0 "Starting mass flowrate";
        parameter Types.Temperature Tstart = 300 "Starting temperature";
        parameter Real sf = 1e2 "Smoothing factor (higher->steeper)";
        constant Real pi = Modelica.Constants.pi;
        Types.MassFlowRate w;
        CombinedCycle.Substances.Gas gas(pstart = pstart, Tstart = Tstart);
        CombinedCycle.Connectors.FlangeB flange(w(start = wstart, nominal = w_min)) annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput load "Percentage of usage" annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        Real lambda(start = 1) "One above load_min_w, zero below load_min_w, with smoothing";
        Modelica.Blocks.Interfaces.RealInput disturbance "Relative change in T and flow rate" annotation(Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = -90, origin = {0, 60}), iconTransformation(extent = {{-20, -20}, {20, 20}}, rotation = -90, origin = {0, 60})));
      equation
        lambda = (atan(sf * (load - load_min_w)) + pi / 2) / pi;
        flange.w = -w;
        w = (w_min * (1 - lambda) + (w_min + (w_max - w_min) * (load - load_min_w) * 2) * lambda) * (1 + disturbance);
        gas.T = (T_min + (T_max - T_min) * (load * 2)) * (1 - lambda) + T_max * lambda + (T_max - T_min) * disturbance;
        flange.p = gas.p;
        flange.h = gas.h;
        annotation(Icon(graphics), Documentation(info = "<html>
<p><b>Modelling options</b></p>
<p>The actual gas used in the component is determined by the replaceable <tt>Medium</tt> package. In the case of multiple component, variable composition gases, the nominal gas composition is given by <tt>Xnom</tt>,whose default value is <tt>Medium.reference_X</tt> .
<p>If <tt>G</tt> is set to zero, the flowrate source is ideal; otherwise, the outgoing flowrate decreases proportionally to the outlet pressure.</p>
<p>If the <tt>in_w0</tt> connector is wired, then the source massflowrate is given by the corresponding signal, otherwise it is fixed to <tt>w0</tt>.</p>
<p>If the <tt>in_T</tt> connector is wired, then the source temperature is given by the corresponding signal, otherwise it is fixed to <tt>T</tt>.</p>
<p>If the <tt>in_X</tt> connector is wired, then the source massfraction is given by the corresponding signal, otherwise it is fixed to <tt>Xnom</tt>.</p>
</html>", revisions = "<html>
<ul>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>w0fix</tt> and <tt>Tfix</tt> and <tt>Xfix</tt>; the connection of external signals is now detected automatically.</li> <br> Adapted to Modelica.Media
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"), Diagram(graphics));
      end BraytonJoule_exhaustGases_disturbances;

      model FlowPConstNvolumes "Heat exchanger for gas including heat transfer model"
        extends CombinedCycle.Icons.Gas.Tube;
        import CombinedCycle.Types.*;
        parameter Integer N(min = 1) = 3 "Number of temperature nodes in the Tube";
        parameter Boolean Direction = true "Used to set parallel or countercurrent heat exchanger";
        parameter Pressure pstart = 101325 "Pressure start value of outgoing gas";
        parameter Temperature Tstart = 300 "Temperature start value of outgoing gas";
        parameter Temperature Tmstart = 300 "Metal wall start temperature";
        parameter Real V = 1 "Inner volume";
        parameter Modelica.SIunits.ThermalConductance Gn = 1 "Thermal conductance between fluid and wall in nominal conditions";
        parameter CombinedCycle.Types.MassFlowRate wn = 1 "Mass flow rate in nominal conditions";
        parameter Real tau = 0.6 "Exponent for the evaluetion of the thermal conductance";
        ThermalConductance G(start = Gn) "Thermal conductance between fluid and wall";
        ThermalConductance Gi(start = Gn / N) "Thermal conductance for each finite volume";
        Temperature Tmean[N](each start = Tmstart) "Mean temperature of the mass flow(each T(start = Tmstart))";
        Real M[N](each start = 30 * V) "Gas total mass";
        CombinedCycle.Substances.Gas gas[N + 1](each pstart = pstart, each Tstart = Tstart);
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Connectors.MultiHeatPort wall(N = N) annotation(Placement(transformation(extent = {{-20, 42}, {20, 82}}, rotation = 0), iconTransformation(extent = {{-40, 40}, {40, 60}})));
      equation
        //Definitions
        G = Gn * (-outlet.w / wn) ^ tau;
        Gi = G / N;
        0 = inlet.w + outlet.w "Mass balance";
        inlet.p = outlet.p "Momentum balance";
        for i in 1:N loop
          M[i] = gas[i + 1].d * V;
          if Direction == true then
            der(gas[i + 1].T) * (gas[i + 1].cv * M[i]) = inlet.w * gas[i].h + outlet.w * gas[i + 1].h + wall.hp[i].Q_flow "Energy balance";
            wall.hp[i].Q_flow = Gi * (wall.hp[i].T - Tmean[i]);
          else
            der(gas[i + 1].T) * (gas[i + 1].cv * M[i]) = inlet.w * gas[i].h + outlet.w * gas[i + 1].h + wall.hp[N - i + 1].Q_flow "Energy balance";
            wall.hp[N - i + 1].Q_flow = Gi * (wall.hp[N - i + 1].T - Tmean[i]);
          end if;
          Tmean[i] = (gas[i].T + gas[i + 1].T) / 2;
          gas[i + 1].p = outlet.p;
        end for;
        // Boundary conditions
        inlet.p = gas[1].p;
        inlet.h = gas[1].h;
        outlet.h = gas[N + 1].h;
        annotation(Placement(transformation(extent = {{-60, 40}, {60, 60}}, rotation = 0)), Icon(graphics = {Rectangle(extent = {{-50, 4}, {32, 0}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, 12}, {32, -8}, {52, 2}, {32, 12}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Documentation(info = "<HTML>
<p>Heat exchanger pipes that use a perfect gas as a fluid.
<p>The heat transfer model assumes that the thermal conductance follows a power law on the mass flow rate.
<p>The parameter N is used to decide how many temperature nodes are used for the heat exchange while the parameter Direction can be used to set up counter current heat exchange.
<p>Pressure loss are considerd negligible, the mass inside each volume is considered constant, the density part of the accumulation term of the energy balance is ignored.
</HTML>"));
      end FlowPConstNvolumes;
    end Gas;

    package Water
      model Turbine
        extends CombinedCycle.Icons.Gas.Turbine;
        parameter Real k = 10;
        parameter CombinedCycle.Types.MassFlowRate w_start = k * 1e5;
        CombinedCycle.Connectors.FlangeA inlet(w(start = w_start)) annotation(Placement(visible = true, transformation(origin = {-68, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-71, 79}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      equation
        inlet.w = k * inlet.p;
        annotation(Diagram(graphics), Icon(graphics), Documentation(info = "<html>
<p>Turbine inlet modeled with a simplified Stodola Ellipse equation.
</html>"));
      end Turbine;

      model TurbineThermalPort
        extends CombinedCycle.Components.Water.Turbine;
        parameter Real Gamma = 1000 "Coefficent of heat flux transfer between vapour and shaft wall";
        parameter CombinedCycle.Types.Temperature Tstart = 800;
        parameter CombinedCycle.Types.Pressure pstart = 3e5;
        parameter CombinedCycle.Types.HeatFlux HF_start = Gamma * (Tstart - 750);
        parameter Types.SpecificEnthalpy hout = 2e6 "Exhaust specific enthalpy";
        Types.HeatFlux HF(start = HF_start);
        CombinedCycle.Substances.WaterVapour wat_vap(Tstart = Tstart, pstart = pstart);
        Connectors.HFP shaft_wall(T(start = Tstart)) annotation(Placement(transformation(extent = {{-10, -52}, {10, -32}}, rotation = 0)));
        CombinedCycle.Blocks.Interfaces.RealOutput TIT annotation(Placement(visible = true, transformation(origin = {4, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-42, 92}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      equation
        //Boundary conditions
        wat_vap.p = inlet.p;
        wat_vap.h = inlet.h;
        HF = Gamma * (wat_vap.T - shaft_wall.T);
        shaft_wall.phi = -HF;
        TIT = wat_vap.T;
        annotation(Diagram(graphics), Icon(graphics), Documentation(info = "<html>
<p>Turbine with added an output for the temperature of the vapour entaering the turbine and a thermal port for the heat exchange with the shaft of the turbine.
<p>The external temperature of the shaft is calculated assuming a fixed hat flux transfer coefficient between the vapour and the metal wall.
</html>"));
      end TurbineThermalPort;

      model SourceW_waterVapour_Input "Flowrate source for gas flows"
        extends CombinedCycle.Icons.Water.SourceW;
        import Modelica.SIunits.*;
        parameter Pressure pstart = 101325 "Start pressure";
        parameter Real Tnom = 300 "Fixed temperature";
        parameter MassFlowRate wstart "Start mass flow rate";
        CombinedCycle.Substances.WaterVapour wat_vap(pstart = pstart, Tstart = Tnom);
        CombinedCycle.Connectors.FlangeB flange(w(start = wstart)) annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput w annotation(Placement(transformation(origin = {-62, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 270)));
      equation
        flange.w = -w;
        wat_vap.T = Tnom "Temperature set by parameter";
        flange.p = wat_vap.p;
        flange.h = wat_vap.h;
        annotation(Icon(graphics), Documentation(info = "<html>
<p>Source of vapour water, the temperature is imposed by a parameter while the massflow is set by an input signal.
</html>"));
      end SourceW_waterVapour_Input;

      model SourceW_waterLiquid_Input "Flowrate source for gas flows"
        extends CombinedCycle.Icons.Water.SourceW;
        parameter Types.Pressure pstart = 101325 "Start pressure";
        parameter Types.Temperature Tnom = 300 "Nominal temperature";
        parameter Types.MassFlowRate wstart = 70 "Start mass flow rate";
        CombinedCycle.Substances.WaterLiquid wat_liq(pstart = pstart, Tstart = Tnom);
        CombinedCycle.Connectors.FlangeB flange(w(start = wstart)) annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput w annotation(Placement(visible = true, transformation(origin = {-62, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 270), iconTransformation(origin = {-60, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 270)));
      equation
        flange.w = -w;
        wat_liq.T = Tnom "Temperature set by parameter";
        flange.p = wat_liq.p;
        flange.h = wat_liq.h;
        annotation(Icon(graphics), Documentation(info = "<html>
<p>Source of liquid water, the temperature is imposed by a parameter while the massflow is set by an input signal.
</html>"));
      end SourceW_waterLiquid_Input;

      model Flow1DVapourPConstNvolumes
        extends CombinedCycle.Icons.Water.Tube;
        import CombinedCycle.Types.*;
        parameter Integer N(min = 1) = 3 "Number of temperature nodes in the Tube";
        parameter Boolean Direction = true "Used to set parallel or countercurrent heat exchanger";
        parameter Pressure pstart = 101325 "Pressure start value of outgoing vapour";
        parameter Temperature Tstart = 300 "Temperature start value of outgoing vapour";
        parameter Temperature Tmstart = 300 "Metal wall start temperature";
        parameter Real V = 1 "Inner volume";
        parameter Real Cm = 1 "Metal Heat Capacity" annotation(Evaluate = true);
        Temperature Tmean[N](each start = Tmstart) "Mean temperature of the mass flow(each T(start = Tmstart))";
        Real M[N](each start = 30 * V) "Vapour total mass in each finite volume";
        CombinedCycle.Substances.WaterVapour vapour[N + 1](each pstart = pstart, each Tstart = Tstart);
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Connectors.MultiHeatPort wall(N = N) annotation(Placement(transformation(extent = {{-20, 42}, {20, 82}}, rotation = 0), iconTransformation(extent = {{-40, 40}, {40, 60}})));
      equation
        //Definitions
        0 = inlet.w + outlet.w "Mass balance";
        inlet.p = outlet.p "Momentum balance";
        for i in 1:N loop
          M[i] = vapour[i + 1].d * V;
          if Direction == true then
            der(vapour[i + 1].T) * Cm / N = inlet.w * vapour[i].h + outlet.w * vapour[i + 1].h + wall.hp[i].Q_flow "Energy balance";
            wall.hp[i].T = Tmean[i];
          else
            der(vapour[i + 1].T) * Cm / N = inlet.w * vapour[i].h + outlet.w * vapour[i + 1].h + wall.hp[N - i + 1].Q_flow "Energy balance";
            wall.hp[N - i + 1].T = Tmean[i];
          end if;
          Tmean[i] = (vapour[i].T + vapour[i + 1].T) / 2;
          vapour[i + 1].p = outlet.p;
        end for;
        // Boundary conditions
        inlet.p = vapour[1].p;
        inlet.h = vapour[1].h;
        outlet.h = vapour[N + 1].h;
        annotation(Placement(transformation(extent = {{-60, 40}, {60, 60}}, rotation = 0)), Icon(graphics = {Rectangle(extent = {{-50, 4}, {32, 0}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, 12}, {32, -8}, {52, 2}, {32, 12}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Documentation(info = "<HTML>
<p>Heat exchanger pipes that use vapour water as a fluid.
<p>The parameter N is used to decide how many temperature nodes are used for the heat exchange while the parameter Direction can be used to set up counter current heat exchange.
<p>Pressure loss are considerd negligible, the thermal inertia of the fluid is considerd negligible compared to that of the metal walls.
</HTML>"));
      end Flow1DVapourPConstNvolumes;

      model Flow1DLiquidPConstNvolumes
        extends CombinedCycle.Icons.Water.Tube;
        import CombinedCycle.Types.*;
        parameter Integer N(min = 1) = 3 "Number of temperature nodes in the Tube";
        parameter Boolean Direction = true "Used to set parallel or countercurrent heat exchanger";
        parameter Pressure pstart = 101325 "Pressure start value of outgoing liquid";
        parameter Temperature Tstart = 300 "Temperature start value of outgoing liquid";
        parameter Temperature Tmstart = 300 "Metal wall start temperature";
        parameter Real V = 1 "Inner volume";
        parameter Real Cm = 1 "Metal Heat Capacity" annotation(Evaluate = true);
        Temperature Tmean[N](each start = Tmstart) "Mean temperature of the mass flow(each T(start = Tmstart))";
        Real M[N](each start = 30 * V) "Liquid total mass in each finite volume";
        CombinedCycle.Substances.WaterLiquid liquid[N + 1](each pstart = pstart, each Tstart = Tstart);
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(transformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(transformation(extent = {{80, -20}, {120, 20}}, rotation = 0)));
        Connectors.MultiHeatPort wall(N = N) annotation(Placement(transformation(extent = {{-20, 42}, {20, 82}}, rotation = 0), iconTransformation(extent = {{-40, 40}, {40, 60}})));
      equation
        //Definitions
        0 = inlet.w + outlet.w "Mass balance";
        inlet.p = outlet.p "Momentum balance";
        for i in 1:N loop
          M[i] = liquid[i + 1].d * V;
          if Direction == true then
            der(liquid[i + 1].T) * Cm / N = inlet.w * liquid[i].h + outlet.w * liquid[i + 1].h + wall.hp[i].Q_flow "Energy balance";
            wall.hp[i].T = Tmean[i];
          else
            der(liquid[i + 1].T) * Cm / N = inlet.w * liquid[i].h + outlet.w * liquid[i + 1].h + wall.hp[N - i + 1].Q_flow "Energy balance";
            wall.hp[N - i + 1].T = Tmean[i];
          end if;
          Tmean[i] = (liquid[i].T + liquid[i + 1].T) / 2;
          liquid[i + 1].p = outlet.p;
        end for;
        // Boundary conditions
        inlet.p = liquid[1].p;
        inlet.h = liquid[1].h;
        outlet.h = liquid[N + 1].h;
        annotation(Placement(transformation(extent = {{-60, 40}, {60, 60}}, rotation = 0)), Icon(graphics = {Rectangle(extent = {{-50, 4}, {32, 0}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{32, 12}, {32, -8}, {52, 2}, {32, 12}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Documentation(info = "<HTML>
<p>Heat exchanger pipes that use liquid water as a fluid.
<p>The parameter N is used to decide how many temperature nodes are used for the heat exchange while the parameter Direction can be used to set up counter current heat exchange.
<p>Pressure loss are considerd negligible, the thermal inertia of the fluid is considerd negligible compared to that of the metal walls.
</HTML>"));
      end Flow1DLiquidPConstNvolumes;

      model DrumNvolumes
        extends CombinedCycle.Icons.Water.Drum;
        import CombinedCycle.Types.*;
        parameter Integer N(min = 1) = 3 "Number of nodes in the heat exchanger";
        parameter Types.Pressure pstart = 101325 "Pressure start value of outgoing gas";
        parameter LiquidMass M_start = Mv_start + Ml_start;
        parameter VapourMass Mv_start = V * alphastart * sat.d_vs_start;
        parameter LiquidMass Ml_start = V * (1 - alphastart) * sat.d_ls_start;
        parameter Real Vv_start = V - Vl_start;
        parameter Real Vl_start = 1;
        parameter Real alphastart = 0.5;
        parameter Real V = 15 "Inner volume";
        parameter Real Cm = 1 "Metal Heat Capacity" annotation(Evaluate = true);
        LiquidMass M(start = M_start) "Total mass";
        VapourMass Mv(start = Mv_start) "Vapour mass";
        LiquidMass Ml(start = Ml_start) "Liquid mass";
        Real Vv(start = Vv_start, nominal = V / 2) "Vapour volume";
        Real Vl(start = Vl_start, nominal = V / 2) "Liquid volume";
        Real alpha(start = alphastart) "Liquid/Vapour volumes ratio";
        Pressure p(start = pstart);
        CombinedCycle.Substances.WaterSaturation sat(pstart = pstart);
        HeatFlowRate Q_tot "Total heat flow through the wall";
        //Dummy variables
        Real der_E(nominal = 1e6);
        Real der_M;
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(transformation(extent = {{-100, -78}, {-60, -38}}, rotation = 0), iconTransformation(extent = {{-100, -78}, {-60, -38}})));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(transformation(extent = {{40, 40}, {80, 80}}, rotation = 0), iconTransformation(extent = {{40, 40}, {80, 80}})));
        Connectors.MultiHeatPort wall(hp(each T(start = sat.T_s_start)), N = N) annotation(Placement(transformation(extent = {{-18, -96}, {22, -56}}, rotation = 0), iconTransformation(extent = {{-40, -100}, {40, -78}})));
        Blocks.Interfaces.RealOutput voidFraction annotation(Placement(transformation(extent = {{66, -70}, {86, -50}}, rotation = 0)));
      equation
        //HP. liquid and vapour are at saturation condition.
        //Definitions
        Q_tot = sum(wall.hp.Q_flow);
        alpha = Vv / V;
        V = Vl + Vv;
        Mv = V * alpha * sat.d_vs;
        Ml = V * (1 - alpha) * sat.d_ls;
        M = Mv + Ml;
        der_M = V * ((-der(alpha) * sat.d_ls) + (1 - alpha) * sat.d_d_ls_dp * der(p) + der(alpha) * sat.d_vs + alpha * sat.d_d_vs_dp * der(p)) "Mass balance";
        der_E = V * ((-der(alpha) * sat.d_ls * sat.u_ls) + (1 - alpha) * sat.d_de_ls_dp * der(p) + der(alpha) * sat.d_vs * sat.u_vs + alpha * sat.d_de_vs_dp * der(p)) + Cm * sat.d_T_s_dp * der(p) "Energy balance";
        der_M = inlet.w + outlet.w "Mass balance";
        der_E = inlet.w * inlet.h + outlet.w * outlet.h + Q_tot "Energy balance";
        inlet.p = outlet.p "Momentum balance";
        inlet.p = p;
        p = sat.p;
        // Boundary conditions
        outlet.h = sat.h_vs;
        wall.hp.T = ones(N) * sat.T_s;
        voidFraction = alpha;
        assert(alpha < 1 and alpha > 0, "The vapour volume shouldn't be more of the total volume or less than zero.");
        annotation(Placement(transformation(extent = {{-60, 80}, {60, 100}}, rotation = 0)), Icon(graphics), Documentation(info = "<html>
<p>This model describes a constant volume buffer with metal walls. The metal wall temperature and the heat transfer coefficient between the wall and the fluid are uniform. The wall is thermally insulated from the outside.</p>
<p>If the inlet or the outlet are connected to a bank of tubes, the model can actually represent a collector or a distributor.</p>
</html>", revisions = "<html>
<ul>
<li><i>30 May 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Initialisation support added.</li>
<li><i>19 Nov 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>5 Mar 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>
"), Diagram(graphics));
      end DrumNvolumes;

      model Flow1DLiquid_inputw
        extends CombinedCycle.Icons.Water.Tube;
        parameter Types.MassFlowRate wstart = 0 "Start mass flow rate";
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(visible = true, transformation(origin = {-100, 1.33227e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -1.33227e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(visible = true, transformation(origin = {100, -1.33227e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {100, -1.33227e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealInput w(start = wstart) annotation(Placement(visible = true, transformation(origin = {-30, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {1.33227e-15, 60}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      equation
        inlet.w = w "control signal";
        inlet.w + outlet.w = 0 "mass balance";
        inlet.h = outlet.h;
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Text(origin = {-1, 87}, extent = {{9, 9}, {-9, -9}}, textString = "w"), Rectangle(origin = {-15, 2}, fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-35, 2}, {35, -2}}), Polygon(origin = {29.72, 2}, fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-9.72361, 10}, {-9.72361, -10}, {10.2764, 0}, {-9.72361, 10}})}), Documentation(info = "<html>
  <p>Fluid valve with a control input. The control signal is the allowed mass flow and not the pressure loss.</p>
  </html>"));
      end Flow1DLiquid_inputw;

      model FlowSplit
        extends CombinedCycle.Icons.Water.FlowSplit;
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -1.33227e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet1 annotation(Placement(visible = true, transformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet2 annotation(Placement(visible = true, transformation(origin = {100, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      equation
        inlet.w + outlet1.w + outlet2.w = 0 "Mass Balance";
        outlet1.p = inlet.p;
        outlet2.p = inlet.p;
        outlet1.h = inlet.h;
        outlet2.h = inlet.h;
        annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Documentation(info = "<html>
  <p>Component to split one inflow in two outflows.</p>
  </html>"));
      end FlowSplit;

      model FlowJoin
        extends CombinedCycle.Icons.Water.FlowJoin;
        CombinedCycle.Connectors.FlangeA inlet1 annotation(Placement(visible = true, transformation(origin = {-50, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeA inlet2 annotation(Placement(visible = true, transformation(origin = {-70, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(visible = true, transformation(origin = {60, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -1.33227e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Substances.WaterVapour wat_vap_out;
        Modelica.Blocks.Interfaces.RealOutput Tout annotation(Placement(visible = true, transformation(origin = {28, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      equation
        inlet1.w + inlet2.w + outlet.w = 0 "Mass Balance";
        inlet1.p = outlet.p;
        inlet2.p = outlet.p;
        outlet.w * outlet.h + inlet1.w * inlet1.h + inlet2.w * inlet2.h = 0 "Energy Balance";
        wat_vap_out.p = outlet.p;
        wat_vap_out.h = outlet.h;
        Tout = wat_vap_out.T;
        annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Documentation(info = "<html>
  <p>Component to calculate the characteristics of a single outflow of vapour water as the mix of two inflows. The Temperature of the outflow is exported as an output</p>
  </html>"));
      end FlowJoin;

      model Valve_DeltaPconst
        extends CombinedCycle.Icons.Water.Tube;
        CombinedCycle.Connectors.FlangeA inlet annotation(Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Connectors.FlangeB outlet annotation(Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        parameter CombinedCycle.Types.Pressure deltaP = 2e5 "Fixed Pressure loss";
      equation
        inlet.w + outlet.w = 0;
        inlet.h = outlet.h;
        outlet.p = inlet.p - deltaP;
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Polygon(origin = {-5, 0}, lineColor = {255, 255, 255}, lineThickness = 1, points = {{-49, 28}, {-49, -28}, {57, 28}, {57, -28}, {-49, 28}})}), Documentation(info = "<html>
  <p>Fluid valve with fixed pressure loss.</p>
  </html>"));
      end Valve_DeltaPconst;
    end Water;

    package Thermal
      model CylinderFourierConstant "Thermal model of a hollow cylinder by Fourier's equation - 1 axial node and Nr radial nodes"
        import Modelica.SIunits.*;
        import CombinedCycle.Choices.CylinderFourier.NodeDistribution;
        extends CombinedCycle.Icons.MetalWall;
        //  replaceable model MaterialModel = MaterialProperties.Metals.StandardSteel constrainedby
        //    MaterialProperties.Interfaces.PartialMaterial "Metal model";
        parameter Integer Nr = 3 "Number of radial nodes";
        /* 
                                                                                                                                                                                                                                                                    parameter NodeDistribution nodeDistribution = NodeDistribution.uniform 
                                                                                                                                                                                                                                                                      "Node distribution";
                                                                                                                                                                                                                                                                    */
        parameter Length rint "Internal radius";
        parameter Length rext "External radius";
        parameter Density rho "Density of the material";
        parameter SpecificHeatCapacity c "Specific heat capacity of material";
        parameter ThermalConductivity lambda "Thermal conductivity of the material";
        parameter Temperature Tstartint = 300 "Temperature start value at rint (first node)" annotation(Dialog(tab = "Initialisation"));
        parameter Temperature Tstartext = 300 "Temperature start value at rext (last node)" annotation(Dialog(tab = "Initialisation"));
        /*  
                                                                                                                                                                                                                                                                    parameter Choices.Init.Options initOpt=Choices.Init.Options.noInit 
                                                                                                                                                                                                                                                                      "Initialisation option" annotation(Dialog(tab = "Initialisation"));
                                                                                                                                                                                                                                                                    */
        /*
                                                                                                                                                                                                                                                                    Length r[Nr](each fixed=false) "Node radii";
                                                                                                                                                                                                                                                                    Length r1_2[Nr-1](each fixed=false) "Slice mean radii";
                                                                                                                                                                                                                                                                    Length r_lin[Nr](each fixed=false) "Linearly distributed radii";
                                                                                                                                                                                                                                                                    */
        final parameter Length r[Nr] = r_lin "Node radii";
        final parameter Length r1_2[Nr - 1] = array((r[i + 1] + r[i]) / 2 for i in 1:Nr - 1) "Slice mean radii";
        final parameter Length r_lin[Nr] = linspace(rint, rext, Nr) "Linearly distributed radii";
        Real A[Nr](each fixed = false);
        Real B[Nr](each fixed = false);
        Real C[Nr](each fixed = false);
        Types.Temperature T[Nr](start = linspace(Tstartint, Tstartext, Nr)) "Nodal temperatures";
        Temperature Tm "Mean temperature";
        // MaterialModel metal[Nr] "Metal properties at the nodes";
        Connectors.HFP internalBoundary annotation(Placement(transformation(extent = {{-20, 20}, {20, 40}}, rotation = 0)));
        Connectors.HFP externalBoundary annotation(Placement(transformation(extent = {{-20, -40}, {20, -20}}, rotation = 0)));
      equation
        // Generation of the temperature node distribution
        /*
                                                                                                                                                                                                                                                                    r_lin = linspace(rint,rext,Nr) "Linearly distributed node radii";
                                                                                                                                                                                                                                                                    for i in 1:Nr loop
                                                                                                                                                                                                                                                                      r[i]= r_lin[i] "Uniform distribution of node radii";
                                                                                                                                                                                                                                                                    end for;
                                                                                                                                                                                                                                                                    */
        /* unsupported in JModelica 1.5
                                                                                                                                                                                                                                                                    for i in 1:Nr loop
                                                                                                                                                                                                                                                                      if nodeDistribution == NodeDistribution.uniform then
                                                                                                                                                                                                                                                                        r[i]= r_lin[i] "Uniform distribution of node radii";
                                                                                                                                                                                                                                                                      elseif nodeDistribution == NodeDistribution.thickInternal then
                                                                                                                                                                                                                                                                        r[i]= rint + 1/(rext-rint)*(r_lin[i]-rint)^2 
                                                                                                                                                                                                                                                                          "Quadratically distributed node radii - thickest at rint";
                                                                                                                                                                                                                                                                      elseif nodeDistribution == NodeDistribution.thickExternal then
                                                                                                                                                                                                                                                                        r[i]= rext - 1/(rext-rint)*(rext-r_lin[i])^2 
                                                                                                                                                                                                                                                                          "Quadratically distributed node radii - thickest at rext";
                                                                                                                                                                                                                                                                      elseif nodeDistribution == NodeDistribution.thickBoth then
                                                                                                                                                                                                                                                                        if r_lin[i] <= (rint+rext)/2 then
                                                                                                                                                                                                                                                                          r[i]=  2/(rext-rint)*(r_lin[i]-rint)^2+rint 
                                                                                                                                                                                                                                                                            "Quadratically distributed node radii - thickest at rint";
                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                          r[i]= -2/(rext-rint)*(r_lin[i]-rext)^2+rext 
                                                                                                                                                                                                                                                                            "Quadratically distributed node radii - thickest at rext";
                                                                                                                                                                                                                                                                        end if;
                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                          r[i] = 0;
                                                                                                                                                                                                                                                                        assert(true,"Unsupported NodeDistribution type");
                                                                                                                                                                                                                                                                      end if;
                                                                                                                                                                                                                                                                    end for;
                                                                                                                                                                                                                                                                    */
        /*
                                                                                                                                                                                                                                                                    for i in 1:Nr-1 loop
                                                                                                                                                                                                                                                                      r1_2[i] = (r[i+1]+r[i])/2;
                                                                                                                                                                                                                                                                    end for;
                                                                                                                                                                                                                                                                    */
        // Spatially discretized coefficients of Fourier's equation
        for i in 2:Nr - 1 loop
          A[i] = r1_2[i - 1] / (r[i] * (r[i] - r[i - 1]) * (r1_2[i] - r1_2[i - 1]));
          C[i] = r1_2[i] / (r[i] * (r[i + 1] - r[i]) * (r1_2[i] - r1_2[i - 1]));
          B[i] = (-A[i]) - C[i];
        end for;
        // Not used by Fourier equations
        A[1] = 0;
        B[1] = 0;
        C[1] = 0;
        A[Nr] = 0;
        B[Nr] = 0;
        C[Nr] = 0;
        // Thermal field
        for i in 2:Nr - 1 loop
          rho * c / lambda * der(T[i]) = A[i] * T[i - 1] + B[i] * T[i] + C[i] * T[i + 1] "Fourier's equation";
        end for;
        // Thermal boundary conditions
        internalBoundary.T = T[1];
        externalBoundary.T = T[Nr];
        internalBoundary.phi = -lambda * (T[2] - T[1]) / (r[2] - r[1]);
        externalBoundary.phi = lambda * (T[Nr] - T[Nr - 1]) / (r[Nr] - r[Nr - 1]);
        // Mean temperature
        Tm = 1 / (rext ^ 2 - rint ^ 2) * sum((T[i] * r[i] + T[i + 1] * r[i + 1]) * (r[i + 1] - r[i]) for i in 1:Nr - 1);
        //  Tm = sum(T)/Nr;
      initial equation
        // Initial conditions
        /*

                                                                                                                                                                                                                                                                    if initOpt == Choices.Init.Options.noInit then
                                                                                                                                                                                                                                                                      // do nothing
                                                                                                                                                                                                                                                                    elseif initOpt == Choices.Init.Options.steadyState then
                                                                                                                                                                                                                                                                      der(T[2:Nr-1]) = zeros(Nr-2);
                                                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                                                      assert(false, "Unsupported initialisation option");
                                                                                                                                                                                                                                                                    end if;
                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                    */
        annotation(Icon(graphics = {Text(extent = {{-94, 52}, {-42, 24}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Forward, textString = "Int"), Text(extent = {{-90, -24}, {-42, -50}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Forward, textString = "Ext"), Text(extent = {{-98, -44}, {102, -72}}, lineColor = {191, 95, 0}, textString = "%name")}));
      end CylinderFourierConstant;

      model CylinderThermalStressConstant "Thermal stress model in a hollow cylinder - constant conductivity"
        extends CombinedCycle.Components.Thermal.CylinderFourierConstant;
        import Modelica.SIunits.*;
        Types.Stress extThermalStress "Thermal stress at the external surface (radial only)";
        Types.Stress intThermalStress "Thermal stress at the internal surface (radial only)";
        parameter Real linearExpansionCoefficient(final unit = "1/K");
        parameter Pressure youngModulus;
        parameter Real poissonRatio;
      equation
        // Thermal stresses
        intThermalStress = linearExpansionCoefficient * youngModulus / (1 - poissonRatio) * (Tm - T[1]);
        extThermalStress = linearExpansionCoefficient * youngModulus / (1 - poissonRatio) * (Tm - T[Nr]);
        annotation(Diagram(graphics), Icon(graphics = {Text(extent = {{-80, 20}, {80, -20}}, lineColor = {255, 255, 255}, textString = "sigma")}), Documentation(info = "<html>
This model extends <tt>CylinderFourier</tt> by adding the computation of the three components of the thermal stress at the internal and external surfaces.
</html>", revisions = "<html>
<ul>
<li><i>10 May 2005</i>
    by <a href=\"mailto:luca.bascetta@polimi.it\">Luca Bascetta</a>:<br>
       First release.</li>
</ul>
</html>"));
      end CylinderThermalStressConstant;

      model Wall_fixedT
        parameter Integer N(min = 1) = 3;
        parameter Types.Temperature T = 300;
        parameter Real G = 176517;
        Connectors.MultiHeatPort wall(N = N) annotation(Placement(visible = true, transformation(origin = {3.55271e-15, 20}, extent = {{-40, -40}, {40, 40}}, rotation = 0), iconTransformation(origin = {3.47378e-15, 20}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
      equation
        for i in 1:N loop
          wall.hp[i].Q_flow = G * (wall.hp[i].T - T);
        end for;
        annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(fillPattern = FillPattern.Solid, extent = {{-80, 20}, {80, -20}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {0, -10}, fillPattern = FillPattern.Solid, extent = {{-80, 30}, {80, -10}})}));
      end Wall_fixedT;
    end Thermal;
  end Components;

  package Icons "Icons for ThermoPower library"
    extends Modelica.Icons.Library;

    package Water "Icons for component using water/steam as working fluid"
      extends Modelica.Icons.Library;

      partial model SourceP
        annotation(Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 0, 0}), Text(extent = {{-20, 34}, {28, -26}}, lineColor = {255, 255, 255}, textString = "P"), Text(extent = {{-100, -78}, {100, -106}}, textString = "%name")}));
      end SourceP;

      partial model SourceW
        annotation(Icon(graphics = {Rectangle(extent = {{-80, 40}, {80, -40}}, lineColor = {0, 0, 0}), Polygon(points = {{-12, -20}, {66, 0}, {-12, 20}, {34, 0}, {-12, -20}}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -52}, {100, -80}}, textString = "%name")}));
      end SourceW;

      partial model Tube
        annotation(Icon(graphics = {Rectangle(extent = {{-80, 40}, {80, -40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder)}), Diagram(graphics));
      end Tube;

      partial model Mixer
        annotation(Icon(graphics = {Ellipse(extent = {{80, 80}, {-80, -80}}, lineColor = {0, 0, 0}), Text(extent = {{-100, -84}, {100, -110}}, textString = "%name")}), Diagram(graphics));
      end Mixer;

      partial model Tank
        annotation(Icon(graphics = {Rectangle(extent = {{-60, 60}, {60, -80}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-54, 60}, {54, 12}}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-54, 12}, {54, -72}}, lineColor = {0, 0, 255})}));
      end Tank;

      partial model Valve
        annotation(Icon(graphics = {Line(points = {{0, 40}, {0, 0}}, color = {0, 0, 0}, thickness = 0.5), Polygon(points = {{-80, 40}, {-80, -40}, {0, 0}, {-80, 40}}, lineColor = {0, 0, 0}, lineThickness = 0.5), Polygon(points = {{80, 40}, {0, 0}, {80, -40}, {80, 40}}, lineColor = {0, 0, 0}, lineThickness = 0.5), Rectangle(extent = {{-20, 60}, {20, 40}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid)}), Diagram(graphics));
      end Valve;

      model FlowJoin
        annotation(Diagram(graphics), Icon(graphics = {Polygon(points = {{-40, 60}, {0, 20}, {40, 20}, {40, -20}, {0, -20}, {-40, -60}, {-40, -20}, {-20, 0}, {-40, 20}, {-40, 60}}, lineColor = {0, 0, 0})}));
      end FlowJoin;

      model FlowSplit
        annotation(Diagram(graphics), Icon(graphics = {Polygon(points = {{40, 60}, {0, 20}, {-40, 20}, {-40, -20}, {0, -20}, {40, -60}, {40, -20}, {22, 0}, {40, 20}, {40, 60}}, lineColor = {0, 0, 0})}));
      end FlowSplit;

      model SensThrough
        annotation(Icon(graphics = {Rectangle(extent = {{-40, -20}, {40, -60}}, lineColor = {0, 0, 0}), Line(points = {{0, 20}, {0, -20}}, color = {0, 0, 0}), Ellipse(extent = {{-40, 100}, {40, 20}}, lineColor = {0, 0, 0}), Line(points = {{40, 60}, {60, 60}}), Text(extent = {{-100, -76}, {100, -100}}, textString = "%name")}));
      end SensThrough;

      model SensP
        annotation(Icon(graphics = {Line(points = {{0, 20}, {0, -20}}, color = {0, 0, 0}), Ellipse(extent = {{-40, 100}, {40, 20}}, lineColor = {0, 0, 0}), Line(points = {{40, 60}, {60, 60}}), Text(extent = {{-100, -52}, {100, -86}}, textString = "%name")}));
      end SensP;

      model Drum
        annotation(Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Polygon(points = {{-60, 0}, {-60, -6}, {-58, -16}, {-52, -30}, {-44, -42}, {-38, -46}, {-32, -50}, {-22, -56}, {-16, -58}, {-8, -60}, {-6, -60}, {0, -60}, {6, -60}, {12, -58}, {22, -56}, {30, -52}, {36, -48}, {42, -42}, {48, -36}, {52, -28}, {58, -18}, {60, -8}, {60, 0}, {-60, 0}}, lineColor = {128, 128, 128}), Polygon(points = {{-60, 0}, {-58, 16}, {-50, 34}, {-36, 48}, {-26, 54}, {-16, 58}, {-6, 60}, {0, 60}, {10, 60}, {20, 56}, {30, 52}, {36, 48}, {46, 40}, {52, 30}, {56, 22}, {58, 14}, {60, 6}, {60, 0}, {-60, 0}}, lineColor = {128, 128, 128}, fillColor = {159, 191, 223}, fillPattern = FillPattern.Solid)}));
      end Drum;

      partial model Pump
        annotation(Icon(graphics = {Polygon(points = {{-40, -24}, {-60, -60}, {60, -60}, {40, -24}, {-40, -24}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-60, 80}, {60, -40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.Sphere), Polygon(points = {{-30, 52}, {-30, -8}, {48, 20}, {-30, 52}}, lineColor = {0, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, fillColor = {255, 255, 255}), Text(extent = {{-100, -64}, {100, -90}}, textString = "%name")}));
      end Pump;

      partial model Accumulator
        annotation(Icon(graphics = {Rectangle(extent = {{-60, 80}, {60, -40}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-60, 100}, {60, 60}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-60, -20}, {60, -60}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-52, 94}, {52, 64}}, lineColor = {0, 0, 191}, pattern = LinePattern.None, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-52, 22}, {52, -40}}, lineColor = {0, 0, 191}, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-52, 80}, {52, 20}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-52, -24}, {52, -54}}, lineColor = {0, 0, 191}, pattern = LinePattern.None, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-4, -58}, {4, -86}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-26, -86}, {26, -94}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid), Text(extent = {{-62, -100}, {64, -122}}, textString = "%name"), Polygon(points = {{-74, 86}, {-60, 72}, {-54, 78}, {-68, 92}, {-74, 86}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid)}), Diagram(graphics));
      end Accumulator;

      partial model PumpMech
        annotation(Icon(graphics = {Rectangle(extent = {{54, 28}, {80, 12}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Polygon(points = {{-40, -24}, {-60, -60}, {60, -60}, {40, -24}, {-40, -24}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 0, 191}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-60, 80}, {60, -40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.Sphere), Polygon(points = {{-30, 52}, {-30, -8}, {48, 20}, {-30, 52}}, lineColor = {0, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.HorizontalCylinder, fillColor = {255, 255, 255}), Text(extent = {{-100, -64}, {100, -90}}, textString = "%name")}));
      end PumpMech;

      partial model PressDrop
        annotation(Icon(graphics = {Rectangle(extent = {{-80, 40}, {80, -40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder), Polygon(points = {{-80, 40}, {-42, 40}, {-20, 12}, {20, 12}, {40, 40}, {80, 40}, {80, -40}, {40, -40}, {20, -12}, {-20, -12}, {-40, -40}, {-80, -40}, {-80, 40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {0, 0, 255})}), Diagram(graphics));
      end PressDrop;

      partial model SteamTurbineUnit
        annotation(Icon(graphics = {Line(points = {{14, 20}, {14, 42}, {38, 42}, {38, 20}}, color = {0, 0, 0}, thickness = 0.5), Rectangle(extent = {{-100, 8}, {100, -8}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Polygon(points = {{-14, 48}, {-14, -48}, {14, -20}, {14, 20}, {-14, 48}}, lineColor = {0, 0, 0}, lineThickness = 0.5, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{38, 20}, {38, -20}, {66, -46}, {66, 48}, {38, 20}}, lineColor = {0, 0, 0}, lineThickness = 0.5, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid), Polygon(points = {{-66, 20}, {-66, -20}, {-40, -44}, {-40, 48}, {-66, 20}}, lineColor = {0, 0, 0}, lineThickness = 0.5, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid), Line(points = {{-100, 70}, {-100, 70}, {-66, 70}, {-66, 20}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{-40, 46}, {-40, 70}, {26, 70}, {26, 42}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{-14, -46}, {-14, -70}, {66, -70}, {66, -46}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{66, -70}, {100, -70}}, color = {0, 0, 255}, thickness = 0.5)}), Diagram(graphics));
      end SteamTurbineUnit;

      partial model Header
        annotation(Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {95, 95, 95}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid), Ellipse(extent = {{70, 70}, {-70, -70}}, lineColor = {95, 95, 95}), Text(extent = {{-100, -84}, {100, -110}}, textString = "%name")}), Diagram(graphics));
      end Header;
    end Water;

    partial model HeatFlow
      annotation(Icon(graphics = {Rectangle(extent = {{-80, 20}, {80, -20}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Forward)}));
    end HeatFlow;

    partial model MetalWall
      annotation(Icon(graphics = {Rectangle(extent = {{-80, 20}, {80, -20}}, lineColor = {0, 0, 0}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid)}));
    end MetalWall;

    package Gas "Icons for component using water/steam as working fluid"
      extends Modelica.Icons.Library;

      partial model SourceP
        annotation(Icon(graphics = {Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {128, 128, 128}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Text(extent = {{-20, 34}, {28, -26}}, lineColor = {255, 255, 255}, textString = "P"), Text(extent = {{-100, -78}, {100, -106}}, textString = "%name")}));
      end SourceP;

      partial model SourceW
        annotation(Icon(graphics = {Rectangle(extent = {{-80, 40}, {80, -40}}, lineColor = {128, 128, 128}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{-12, -20}, {66, 0}, {-12, 20}, {34, 0}, {-12, -20}}, lineColor = {128, 128, 128}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -52}, {100, -80}}, textString = "%name")}));
      end SourceW;

      partial model Tube
        annotation(Icon(graphics = {Rectangle(extent = {{-80, 40}, {80, -40}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {159, 159, 223})}), Diagram(graphics));
      end Tube;

      partial model Mixer
        annotation(Icon(graphics = {Ellipse(extent = {{80, 80}, {-80, -80}}, lineColor = {128, 128, 128}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Text(extent = {{-100, -84}, {100, -110}}, textString = "%name")}), Diagram(graphics));
      end Mixer;

      partial model Valve
        annotation(Icon(graphics = {Line(points = {{0, 40}, {0, 0}}, color = {0, 0, 0}, thickness = 0.5), Polygon(points = {{-80, 40}, {-80, -40}, {0, 0}, {-80, 40}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{80, 40}, {0, 0}, {80, -40}, {80, 40}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-20, 60}, {20, 40}}, lineColor = {0, 0, 0}, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid)}), Diagram(graphics));
      end Valve;

      model FlowJoin
        annotation(Diagram(graphics), Icon(graphics = {Polygon(points = {{-40, 60}, {0, 20}, {40, 20}, {40, -20}, {0, -20}, {-40, -60}, {-40, -20}, {-20, 0}, {-40, 20}, {-40, 60}}, lineColor = {128, 128, 128}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid)}));
      end FlowJoin;

      model FlowSplit
        annotation(Diagram(graphics), Icon(graphics = {Polygon(points = {{40, 60}, {0, 20}, {-40, 20}, {-40, -20}, {0, -20}, {40, -60}, {40, -20}, {22, 0}, {40, 20}, {40, 60}}, lineColor = {128, 128, 128}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid)}));
      end FlowSplit;

      model SensThrough
        annotation(Icon(graphics = {Rectangle(extent = {{-40, -20}, {40, -60}}, lineColor = {128, 128, 128}, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Line(points = {{0, 20}, {0, -20}}, color = {0, 0, 0}), Ellipse(extent = {{-40, 100}, {40, 20}}, lineColor = {0, 0, 0}), Line(points = {{40, 60}, {60, 60}}), Text(extent = {{-100, -76}, {100, -100}}, textString = "%name")}));
      end SensThrough;

      model SensP
        annotation(Icon(graphics = {Line(points = {{0, 20}, {0, -20}}, color = {0, 0, 0}), Ellipse(extent = {{-40, 100}, {40, 20}}, lineColor = {0, 0, 0}), Line(points = {{40, 60}, {60, 60}}), Text(extent = {{-100, -52}, {100, -86}}, textString = "%name")}));
      end SensP;

      partial model Compressor
        annotation(Icon(graphics = {Polygon(points = {{24, 26}, {30, 26}, {30, 76}, {60, 76}, {60, 82}, {24, 82}, {24, 26}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{-30, 76}, {-30, 56}, {-24, 56}, {-24, 82}, {-60, 82}, {-60, 76}, {-30, 76}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-60, 8}, {60, -8}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Polygon(points = {{-30, 60}, {-30, -60}, {30, -26}, {30, 26}, {-30, 60}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid)}), Diagram(graphics));
      end Compressor;

      partial model Turbine
        annotation(Icon(graphics = {Polygon(points = {{-28, 76}, {-28, 28}, {-22, 28}, {-22, 82}, {-60, 82}, {-60, 76}, {-28, 76}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{26, 56}, {32, 56}, {32, 76}, {60, 76}, {60, 82}, {26, 82}, {26, 56}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-60, 8}, {60, -8}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Polygon(points = {{-28, 28}, {-28, -26}, {32, -60}, {32, 60}, {-28, 28}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid)}), Diagram(graphics));
      end Turbine;

      partial model GasTurbineUnit
        annotation(Icon(graphics = {Line(points = {{-22, 26}, {-22, 48}, {22, 48}, {22, 28}}, color = {0, 0, 0}, thickness = 2.5), Rectangle(extent = {{-100, 8}, {100, -8}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, fillColor = {160, 160, 164}), Polygon(points = {{-80, 60}, {-80, -60}, {-20, -26}, {-20, 26}, {-80, 60}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Polygon(points = {{20, 28}, {20, -26}, {80, -60}, {80, 60}, {20, 28}}, lineColor = {128, 128, 128}, lineThickness = 0.5, fillColor = {159, 159, 223}, fillPattern = FillPattern.Solid), Ellipse(extent = {{-16, 64}, {16, 32}}, lineColor = {0, 0, 0}, fillPattern = FillPattern.Sphere, fillColor = {255, 0, 0})}), Diagram(graphics));
      end GasTurbineUnit;
    end Gas;
  end Icons;

  package Choices "Choice enumerations for ThermoPower models"
    package CylinderFourier
      type NodeDistribution = enumeration(uniform "Uniform distribution of node radii", thickInternal "Quadratically distributed node radii - thickest at rint", thickExternal "Quadratically distributed node radii - thickest at rext", thickBoth "Quadratically distributed node radii - thickest at both extremes") "Type, constants and menu choices for node distribution";
    end CylinderFourier;

    package CylinderMechanicalStress
      type MechanicalStandard = enumeration(TRDstandard "TRD standard", ASMEstandard "Laborelec-ASME standard") "Types, constants and menu choices for mechanical standard";
    end CylinderMechanicalStress;

    package Flow1D
      type FFtypes = enumeration(Kfnom "Kfnom friction factor", OpPoint "Friction factor defined by operating point", Cfnom "Cfnom friction factor", Colebrook "Colebrook's equation", NoFriction "No friction") "Type, constants and menu choices to select the friction factor";
      type HCtypes = enumeration(Middle "Middle of the pipe", Upstream "At the inlet", Downstream "At the outlet") "Type, constants and menu choices to select the location of the hydraulic capacitance";
    end Flow1D;

    package PressDrop
      type FFtypes = enumeration(Kf "Kf friction factor", OpPoint "Friction factor defined by operating point", Kinetic "Kinetic friction factor") "Type, constants and menu choices to select the friction factor";
    end PressDrop;

    package Valve
      type CvTypes = enumeration(Av "Av (metric) flow coefficient", Kv "Kv (metric) flow coefficient", Cv "Cv (US) flow coefficient", OpPoint "Av defined by nominal operating point") "Type, constants and menu choices to select the type of Cv data";
    end Valve;

    package TurboMachinery
      type TableTypes = enumeration(matrix "Explicitly supplied as parameter matrix table", file "Read from a file") "Type, constants and menu choices to select the representation of table matrix";
    end TurboMachinery;

    package Init "Options for initialisation"
      type Options = enumeration(noInit "No initial equations", steadyState "Steady-state initialisation", steadyStateNoP "Steady-state initialisation except pressures", steadyStateNoT "Steady-state initialisation except temperatures", steadyStateNoPT "Steady-state initialisation except pressures and temperatures") "Type, constants and menu choices to select the initialisation options";
    end Init;

    package FlowReversal "Options for flow reversal support"
      type Options = enumeration(fullFlowReversal "Full flow reversal support", smallFlowReversal "Small flow reversal allowed (approx. model)", noFlowReversal "Flow reversal is not allowed") "Type, constants and menu choices to select the flow reversal support options";
    end FlowReversal;

    package System
      type Dynamics = enumeration(DynamicFreeInitial "DynamicFreeInitial -- Dynamic balance, Initial guess value", FixedInitial "FixedInitial -- Dynamic balance, Initial value fixed", SteadyStateInitial "SteadyStateInitial -- Dynamic balance, Steady state initial with guess value", SteadyState "SteadyState -- Steady state balance, Initial guess value") "Enumeration to define definition of balance equations";
    end System;

    package FluidPhase
      type FluidPhases = enumeration(Liquid "Liquid", Steam "Steam", TwoPhases "Two Phases") "Type, constants and menu choices to select the fluid phase";
    end FluidPhase;
  end Choices;

  package Functions "Miscellaneous functions"
    extends Modelica.Icons.Library;

    function linear
      extends Modelica.Icons.Function;
      input Real x;
      output Real y;
    algorithm
      y := x;
      annotation(derivative = Functions.linear_der);
    end linear;

    function linear_der
      extends Modelica.Icons.Function;
      input Real x;
      input Real der_x;
      output Real der_y;
    algorithm
      der_y := der_x;
    end linear_der;

    function one
      extends Modelica.Icons.Function;
      input Real x;
      output Real y;
    algorithm
      y := 1;
      annotation(derivative = Functions.one_der);
    end one;

    function one_der
      extends Modelica.Icons.Function;
      input Real x;
      input Real der_x;
      output Real der_y;
    algorithm
      der_y := 0;
    end one_der;

    function sqrtReg "Symmetric square root approximation with finite derivative in zero"
      extends Modelica.Icons.Function;
      input Real x;
      input Real delta = 0.01 "Range of significant deviation from sqrt(x)";
      output Real y;
    algorithm
      y := x / sqrt(sqrt(x * x + delta * delta));
      annotation(derivative(zeroDerivative = delta) = Functions.sqrtReg_der, Documentation(info = "<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/sqrt(delta)</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 16% around x=0.1, 0.25% around x=0.1 and 0.0025% around x=1.
</p> 
</html>", revisions = "<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"), Documentation(info = "<html>
This function approximates sqrt(x)*sign(x), such that the derivative is finite and smooth in x=0. 
</p>
<p>
<table border=1 cellspacing=0 cellpadding=2> 
<tr><th>Function</th><th>Approximation</th><th>Range</th></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= sqrt(abs(x))*sign(x)</td><td>abs(x) &gt;&gt delta</td></tr>
<tr><td>y = sqrtReg(x)</td><td>y ~= x/delta</td><td>abs(x) &lt;&lt  delta</td></tr>
</table>
<p>
With the default value of delta=0.01, the difference between sqrt(x) and sqrtReg(x) is 0.5% around x=0.1 and 0.005% around x=1.
</p> 
</html>", revisions = "<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
    end sqrtReg;

    function sqrtReg_der "Derivative of sqrtReg"
      extends Modelica.Icons.Function;
      input Real x;
      input Real delta = 0.01 "Range of significant deviation from sqrt(x)";
      input Real dx "Derivative of x";
      output Real dy;
    algorithm
      dy := dx * 0.5 * (x * x + 2 * delta * delta) / (x * x + delta * delta) ^ 1.25;
      annotation(Documentation(info = "<html>
</html>", revisions = "<html>
<ul>
<li><i>15 Mar 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Created. </li>
</ul>
</html>"));
    end sqrtReg_der;

    block OffsetController "Offset computation for steady-state conditions"
      extends Modelica.Blocks.Interfaces.BlockIcon;
      parameter Real steadyStateGain = 0.0 "0.0: Adds offset to input - 1.0: Closed loop action to find steady state";
      parameter Real SP0 "Initial setpoint for the controlled variable";
      parameter Real deltaSP = 0 "Variation of the setpoint for the controlled variable";
      parameter Modelica.SIunits.Time Tstart = 0 "Start time of the setpoint ramp change";
      parameter Modelica.SIunits.Time Tend = 0 "End time of the setpoint ramp change";
      parameter Real Kp "Proportional gain";
      parameter Real Ti "Integral time constant";
      parameter Real biasCO "Bias value of the control variable when computing the steady state";
    protected
      Real SP;
      Real error;
      Real integralError;
    public
      Modelica.Blocks.Interfaces.RealInput deltaCO annotation(Placement(transformation(extent = {{-140, 62}, {-100, 100}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput PV annotation(Placement(transformation(extent = {{-140, -100}, {-100, -60}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput CO annotation(Placement(transformation(extent = {{100, -20}, {140, 20}}, rotation = 0)));
    equation
      SP = if time <= Tstart then SP0 else if time >= Tend then SP0 + deltaSP else SP0 + (time - Tstart) / (Tend - Tstart) * deltaSP;
      error = (SP - PV) * steadyStateGain;
      der(integralError) = error;
      CO = Kp * (error + integralError / Ti) + biasCO + (1.0 - steadyStateGain) * deltaCO;
      annotation(Documentation(info = "<HTML>
<p>This model is useful to compute the steady state value of a control variable corresponding to some specified setpoint of an output variable, and to reuse it later to perform simulations starting from this steady state condition.
<p>The block has two different behaviours, depending on the value of the <tt>steadyState</tt> parameter.
<p>When <tt>steadyState = 1</tt>, the <tt>deltaCO</tt> input is ignored, and the block acts as a standard PI controller with transfer function Kp*(1+1/sTi) to bring the process variable connected to the <tt>PV</tt> input at the setpoint value, by acting on the control variable connected to the <tt>CO</tt> output. The setpoint value is <tt>SP0</tt> at time zero, and may change by <tt>deltaSP</tt> from <tt>Tstart</tt> to <tt>Tend</tt>; this can be useful to bring the process far away from the tentative start values of the transient without any inconvenience. The control variable can be biased by <tt>biasCO</tt> to start near the expected steady state value of <tt>CO</tt>.
<p>When <tt>steadyState = 0</tt>, the <tt>PV</tt> input is ignored, and the <tt>CO</tt> output is simply the sum of the <tt>deltaCO</tt> input and of the freezed steady-state output of the controller.
<p>To perform a steady state computation:
<ol>
<li>Set <tt>steadyState = 1</tt> and suitably tune <tt>Kp</tt>, <tt>Ti</tt> and <tt>biasCO</tt>
<li>Simulate a transient until the desired steady state is achieved.
<li>Set <tt>steadyState = 0</tt> and continue the simulation for 0 s
<li>Save the final state of the simulation, which contains the initial steady-state values of all the variables for subsequent transient simulations
</ol>
<p>To perform experiments starting from a steady state:
<ol>
<li>Load a previously saved steady state, to be used as initial state
<li>Perform the simulation of the desired transient. The <tt>offsetCO</tt> input value will be automatically added to the previously computed steady state value.
</ol>
<p><b>Revision history:</b></p>
<ul>
<li><i>15 Feb 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</HTML>"), Diagram(graphics), Icon(graphics = {Text(extent = {{-90, 90}, {94, -92}}, textString = "SS Offset")}));
    end OffsetController;
    annotation(Documentation(info = "<HTML>
This package contains general-purpose functions and models
</HTML>"));
  end Functions;

  package Utils
    model VapourCorrelation
      import Modelica.SIunits.*;
      import Poly = Modelica.Media.Incompressible.TableBased.Polynomials_Temp;
      import Modelica.Utilities.Streams.*;
      parameter Pressure p_min = 20e5;
      parameter Pressure p_max = 100e5;
      parameter Pressure p_nom = (p_max + p_min) / 2;
      parameter Integer N = 50;
      parameter Integer ord_T = 3;
      parameter Integer ord_rhols = 3;
      parameter Integer ord_rhovs = 3;
      parameter Integer ord_rels = 3;
      parameter Integer ord_revs = 3;
      parameter Integer ord_hls = 3;
      parameter Integer ord_hvs = 3;
      parameter Integer ord_drhols_dp = ord_rhols - 1;
      parameter Integer ord_drhovs_dp = ord_rhovs - 1;
      parameter Integer ord_drels_dp = ord_rels - 1;
      parameter Integer ord_drevs_dp = ord_revs - 1;
      package Medium = Modelica.Media.Water.WaterIF97_ph;
      Medium.SaturationProperties sat;
      Pressure p;
      Temperature T;
      Temperature T_approx;
      Density rhols;
      Density rhols_approx;
      Density rhovs;
      Density rhovs_approx;
      Real rels;
      Real rels_approx;
      Real revs;
      Real revs_approx;
      Real hls;
      Real hls_approx;
      Real hvs;
      Real hvs_approx;
      Real dT_dp;
      Real dT_dp_approx;
      Real drhols_dp;
      Real drhols_dp_approx;
      Real drhovs_dp;
      Real drhovs_dp_approx;
      Real drels_dp;
      Real drels_dp_approx;
      Real drevs_dp;
      Real drevs_dp_approx;
      Real r;
      Temperature T_s;
      Real d_ls;
      Real d_vs;
      Real de_ls;
      Real de_vs;
      Real h_ls;
      Real h_vs;
      Real d_T_s_dp;
      Real d_d_ls_dp;
      Real d_d_vs_dp;
      Real d_de_ls_dp;
      Real d_de_vs_dp;
      //Real dummy_rhols;
      //Real dummy_drhols_dp;
    protected
      parameter Pressure p_data[:] = linspace(p_min, p_max, N);
      parameter Density T_data[:] = array(Medium.saturationTemperature(p_data[i]) for i in 1:N);
      parameter Density rhols_data[:] = array(Medium.bubbleDensity(Medium.setSat_p(p_data[i])) for i in 1:N);
      parameter Density rhovs_data[:] = array(Medium.dewDensity(Medium.setSat_p(p_data[i])) for i in 1:N);
      parameter Density rels_data[:] = array(Medium.bubbleDensity(Medium.setSat_p(p_data[i])) * Medium.bubbleEnthalpy(Medium.setSat_p(p_data[i])) - p_data[i] for i in 1:N);
      parameter Density revs_data[:] = array(Medium.dewDensity(Medium.setSat_p(p_data[i])) * Medium.dewEnthalpy(Medium.setSat_p(p_data[i])) - p_data[i] for i in 1:N);
      parameter SpecificEnthalpy hls_data[:] = array(Medium.bubbleEnthalpy(Medium.setSat_p(p_data[i])) for i in 1:N);
      parameter SpecificEnthalpy hvs_data[:] = array(Medium.dewEnthalpy(Medium.setSat_p(p_data[i])) for i in 1:N);
      parameter Real coeff_T[:] = Poly.fitting(p_data / p_nom, T_data, ord_T);
      parameter Real coeff_rhols[:] = Poly.fitting(p_data / p_nom, rhols_data, ord_rhols);
      parameter Real coeff_rhovs[:] = Poly.fitting(p_data / p_nom, rhovs_data, ord_rhovs);
      parameter Real coeff_rels[:] = Poly.fitting(p_data / p_nom, rels_data, ord_rels);
      parameter Real coeff_revs[:] = Poly.fitting(p_data / p_nom, revs_data, ord_revs);
      parameter Real coeff_hls[:] = Poly.fitting(p_data / p_nom, hls_data, ord_hls);
      parameter Real coeff_hvs[:] = Poly.fitting(p_data / p_nom, hvs_data, ord_hvs);
      parameter Real coeff_T_dp[:] = array(coeff_T[i] * (size(coeff_T, 1) - i) for i in 1:size(coeff_T, 1) - 1) / p_nom;
      parameter Real coeff_drhols_dp[:] = array(coeff_rhols[i] * (size(coeff_rhols, 1) - i) for i in 1:size(coeff_rhols, 1) - 1) / p_nom;
      parameter Real coeff_drhovs_dp[:] = array(coeff_rhovs[i] * (size(coeff_rhovs, 1) - i) for i in 1:size(coeff_rhols, 1) - 1) / p_nom;
      parameter Real coeff_drels_dp[:] = array(coeff_rels[i] * (size(coeff_rels, 1) - i) for i in 1:size(coeff_rels, 1) - 1) / p_nom;
      parameter Real coeff_drevs_dp[:] = array(coeff_revs[i] * (size(coeff_revs, 1) - i) for i in 1:size(coeff_revs, 1) - 1) / p_nom;
      //Dummy variables
    equation
      p = p_min + (p_max - p_min) * time;
      sat.psat = p;
      sat.Tsat = Medium.saturationTemperature(p);
      T = sat.Tsat;
      rhols = Medium.bubbleDensity(sat);
      rhovs = Medium.dewDensity(sat);
      rels = Medium.bubbleDensity(sat) * Medium.bubbleEnthalpy(sat) - p;
      revs = Medium.dewDensity(sat) * Medium.dewEnthalpy(sat) - p;
      hvs = Medium.dewEnthalpy(sat);
      hls = Medium.bubbleEnthalpy(sat);
      dT_dp = Medium.saturationTemperature_derp(p);
      drhols_dp = Medium.dBubbleDensity_dPressure(sat);
      drhovs_dp = Medium.dDewDensity_dPressure(sat);
      drels_dp = Medium.dBubbleDensity_dPressure(sat) * Medium.bubbleEnthalpy(sat) + Medium.bubbleDensity(sat) * Medium.dBubbleEnthalpy_dPressure(sat) - 1;
      drevs_dp = Medium.dDewDensity_dPressure(sat) * Medium.dewEnthalpy(sat) + Medium.dewDensity(sat) * Medium.dDewEnthalpy_dPressure(sat) - 1;
      r = hvs - hls;
      T_approx = Poly.evaluate(coeff_T, p / p_nom);
      rhols_approx = Poly.evaluate(coeff_rhols, p / p_nom);
      rhovs_approx = Poly.evaluate(coeff_rhovs, p / p_nom);
      rels_approx = Poly.evaluate(coeff_rels, p / p_nom);
      revs_approx = Poly.evaluate(coeff_revs, p / p_nom);
      hls_approx = Poly.evaluate(coeff_hls, p / p_nom);
      hvs_approx = Poly.evaluate(coeff_hvs, p / p_nom);
      dT_dp_approx = Poly.evaluate(coeff_T_dp, p / p_nom);
      drhols_dp_approx = Poly.evaluate(coeff_drhols_dp, p / p_nom);
      drhovs_dp_approx = Poly.evaluate(coeff_drhovs_dp, p / p_nom);
      drels_dp_approx = Poly.evaluate(coeff_drels_dp, p / p_nom);
      drevs_dp_approx = Poly.evaluate(coeff_drevs_dp, p / p_nom);
      //dummy variables for checking
      T_s = (((+20.31) * p / 6e+006 + (-91.2087)) * p / 6e+006 + 186.107) * p / 6e+006 + 433.8;
      d_ls = (((+(-19.8527)) * p / 6e+006 + 83.716) * p / 6e+006 + (-219.824)) * p / 6e+006 + 913.719;
      d_vs = (((+1.30191) * p / 6e+006 + 0.432043) * p / 6e+006 + 28.7046) * p / 6e+006 + 0.377505;
      de_ls = (((+6.7412e+007) * p / 6e+006 + (-3.09753e+008)) * p / 6e+006 + 5.27973e+008) * p / 6e+006 + 6.29324e+008;
      de_vs = (((+2.38902e+006) * p / 6e+006 + 1.38667e+006) * p / 6e+006 + 7.5263e+007) * p / 6e+006 + 772937;
      h_ls = (((+87156.9) * p / 6e+006 + (-382184)) * p / 6e+006 + 837267) * p / 6e+006 + 672632;
      h_vs = (((+25274.6) * p / 6e+006 + (-125274)) * p / 6e+006 + 108335) * p / 6e+006 + 2.77663e+006;
      d_T_s_dp = ((+1.0155e-005) * p / 6e+006 + (-3.04029e-005)) * p / 6e+006 + 3.10178e-005;
      d_d_ls_dp = ((+(-9.92635e-006)) * p / 6e+006 + 2.79053e-005) * p / 6e+006 + (-3.66374e-005);
      d_d_vs_dp = ((+6.50955e-007) * p / 6e+006 + 1.44014e-007) * p / 6e+006 + 4.7841e-006;
      d_de_ls_dp = ((+33.706) * p / 6e+006 + (-103.251)) * p / 6e+006 + 87.9954;
      d_de_vs_dp = ((+1.19451) * p / 6e+006 + 0.462225) * p / 6e+006 + 12.5438;
      when terminal() then
        //Polynomials
        print("\n\n\nPolynomials:\n");
        //Saturation temperature
        print(PolynomialFunction("T_s", coeff_T, "p/" + String(p_nom)));
        //Liquid density at saturation
        print(PolynomialFunction("d_ls", coeff_rhols, "p/" + String(p_nom)));
        //Vapour density at saturation
        print(PolynomialFunction("d_vs", coeff_rhovs, "p/" + String(p_nom)));
        //Liquid density - specific energy product at saturation
        print(PolynomialFunction("de_ls", coeff_rels, "p/" + String(p_nom)));
        //Vapour density - specific energy product at saturation
        print(PolynomialFunction("de_vs", coeff_revs, "p/" + String(p_nom)));
        //Liquid enthalpy at saturation
        print(PolynomialFunction("h_ls", coeff_hls, "p/" + String(p_nom)));
        //Vapour enthalpy at saturation
        print(PolynomialFunction("h_vs", coeff_hvs, "p/" + String(p_nom)));
        //Derivative of the saturation temperature with respect to the pressure
        print(PolynomialFunction("d_T_s_dp", coeff_T_dp, "p/" + String(p_nom)));
        //Derivative of the liquid density with respect to the pressure
        print(PolynomialFunction("d_d_ls_dp", coeff_drhols_dp, "p/" + String(p_nom)));
        //Derivative of the vapour density with respect to the pressure
        print(PolynomialFunction("d_d_vs_dp", coeff_drhovs_dp, "p/" + String(p_nom)));
        //Derivative of the liquid density - specific energy product at saturation with respect to the pressure
        print(PolynomialFunction("d_de_ls_dp", coeff_drels_dp, "p/" + String(p_nom)));
        //Derivative of the vapour density - specific energy product at saturation with respect to the pressure
        print(PolynomialFunction("d_de_vs_dp", coeff_drevs_dp, "p/" + String(p_nom)));
      end when;
      annotation(uses(Modelica(version = "2.2.1")), experiment, experimentSetupOutput);
    end VapourCorrelation;

    function PolynomialFunction
      input String variable;
      input Real coeff[:];
      input String independent_var;
      output String res;
    protected
      Integer coeff_size;
    algorithm
      coeff_size := size(coeff, 1);
      res := variable + " = ";
      for i in 1:coeff_size - 1 loop
        res := res + "(";
      end for;
      for i in 1:coeff_size loop
        if i < coeff_size then
          res := res + "+(" + String(coeff[i]) + "))*" + independent_var;
        else
          res := res + "+(" + String(coeff[i]) + ");";
        end if;
      end for;
      //    res:= res+"\n";
    end PolynomialFunction;
  end Utils;

  package Blocks
    package Continuous
      block PI "Proportional-Integral controller"
        parameter Real k = 1 "Gain";
        parameter Modelica.SIunits.Time T = 1 "Time Constant (T>0 required)";
        parameter Real x_start = 0 "Initial or guess value of state" annotation(Dialog(group = "Initialization"));
        parameter Real y_start = 0 "Initial value of output" annotation(Dialog(enable = initType == Init.SteadyState or initType == Init.InitialOutput, group = "Initialization"));
        extends CombinedCycle.Blocks.Interfaces.SISO;
        output Real x(start = x_start) "State of block";
      equation
        der(x) = u / T;
        y = k * (x + u);
        annotation(defaultComponentName = "PI", Window(x = 0.24, y = 0.13, width = 0.56, height = 0.5), Documentation(info = "
<HTML>
<p>
This blocks defines the transfer function between the input u and
the output y (element-wise) as <i>PI</i> system:
</p>
<pre>
                 1
   y = k * (1 + ---) * u
                T*s
           T*s + 1
     = k * ------- * u
             T*s
</pre>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general model class <b>TransferFunction</b>
instead and model a PI SISO system with parameters<br>
b = {k*T, k}, a = {T, 0}.
</p>
<pre>
Example:
 
   parameter: k = 0.3,  T = 0.4
 
   results in:
               0.4 s + 1
      y = 0.3 ----------- * u
                 0.4 s
</pre>

<p>
It might be difficult to initialize the PI component in steady state
due to the integrator part.
This is discussed in the description of package
<a href=\"Modelica://Modelica.Blocks.Continuous#info\">Continuous</a>.
</p>
 
</HTML>
"), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-60, 60}, {60, -60}}), Text(extent = {{-68, 24}, {-24, -18}}, lineColor = {0, 0, 0}, textString = "k"), Text(extent = {{-32, 48}, {60, 0}}, lineColor = {0, 0, 0}, textString = "T s + 1"), Text(extent = {{-30, -8}, {52, -40}}, lineColor = {0, 0, 0}, textString = "T s"), Line(points = {{-24, 0}, {54, 0}}, color = {0, 0, 0}), Line(points = {{-100, 0}, {-60, 0}}), Line(points = {{62, 0}, {100, 0}})}), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Line(points = {{-80, 78}, {-80, -90}}, color = {192, 192, 192}), Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 88}, {-80, 90}}), Line(points = {{-90, -80}, {82, -80}}, color = {192, 192, 192}), Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{90, -80}, {68, -72}, {68, -88}, {90, -80}}), Line(points = {{-80, -80}, {-80, -20}, {60, 80}}, color = {0, 0, 255}), Text(lineColor = {192, 192, 192}, extent = {{0, 6}, {60, -56}}, textString = "PI"), Text(extent = {{-150, -150}, {150, -110}}, textString = "T=%T")}));
      end PI;

      block FirstOrder "First order transfer function block (= 1 pole)"
        parameter Real k = 1 "Gain";
        parameter Modelica.SIunits.Time T = 1 "Time Constant";
        parameter Real y_start = 0 "Initial or guess value of output (= state)" annotation(Dialog(group = "Initialization"));
        extends CombinedCycle.Blocks.Interfaces.SISO(y(start = y_start));
      initial equation
        der(y) = 0;
      equation
        der(y) = (k * u - y) / T;
        annotation(Documentation(info = "<HTML>
<p>
This blocks defines the transfer function between the input u
and the output y (element-wise) as <i>first order</i> system:
</p>
<pre>
               k
     y = ------------ * u
            T * s + 1
</pre>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general block <b>TransferFunction</b> instead
and model a first order SISO system with parameters<br>
b = {k}, a = {T, 1}.
</p>
<pre>
Example:
   parameter: k = 0.3, T = 0.4
   results in:
             0.3
      y = ----------- * u
          0.4 s + 1.0
</pre>

</HTML>
"), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Line(points = {{-80, 78}, {-80, -90}}, color = {192, 192, 192}), Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 88}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-90, -80}, {82, -80}}, color = {192, 192, 192}), Polygon(points = {{90, -80}, {68, -72}, {68, -88}, {90, -80}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, -80}, {-70, -45.11}, {-60, -19.58}, {-50, -0.9087}, {-40, 12.75}, {-30, 22.75}, {-20, 30.06}, {-10, 35.41}, {0, 39.33}, {10, 42.19}, {20, 44.29}, {30, 45.82}, {40, 46.94}, {50, 47.76}, {60, 48.36}, {70, 48.8}, {80, 49.12}}), Text(extent = {{0, 0}, {60, -60}}, lineColor = {192, 192, 192}, textString = "PT1"), Text(extent = {{-150, -150}, {150, -110}}, lineColor = {0, 0, 0}, textString = "T=%T")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Text(extent = {{-48, 52}, {50, 8}}, lineColor = {0, 0, 0}, textString = "k"), Text(extent = {{-54, -6}, {56, -56}}, lineColor = {0, 0, 0}, textString = "T s + 1"), Line(points = {{-50, 0}, {50, 0}}, color = {0, 0, 0}), Rectangle(extent = {{-60, 60}, {60, -60}}), Line(points = {{-100, 0}, {-60, 0}}), Line(points = {{60, 0}, {100, 0}})}), Window(x = 0.15, y = 0.04, width = 0.52, height = 0.55));
      end FirstOrder;

      model Actuator
        CombinedCycle.Blocks.Interfaces.RealOutput y annotation(Placement(visible = true, transformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 1.11022e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        parameter Real max(max = 1) = 1;
      equation
        if u < 0 then
          y = 0;
        elseif u > max then
          y = max;
        else
          y = u;
        end if;
        annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {-18, -20}, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-42, -40}, {78, 80}})}));
      end Actuator;

      model I
        parameter Real k = 1 "Gain";
        parameter Real x_start = y_start / k "Initial or guess value of state" annotation(Dialog(group = "Initialization"));
        parameter Real y_start = 0 "Initial value of output" annotation(Dialog(enable = initType == Init.SteadyState or initType == Init.InitialOutput, group = "Initialization"));
        extends CombinedCycle.Blocks.Interfaces.SISO;
        output Real x(start = x_start) "State of block";
      equation
        der(x) = u;
        y = k * x;
        annotation(defaultComponentName = "I", Window(x = 0.24, y = 0.13, width = 0.56, height = 0.5), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Rectangle(extent = {{-60, 60}, {60, -60}}), Text(extent = {{-68, 24}, {-24, -18}}, lineColor = {0, 0, 0}, textString = "k"), Text(extent = {{-32, 48}, {60, 0}}, lineColor = {0, 0, 0}, textString = "T s + 1"), Text(extent = {{-30, -8}, {52, -40}}, lineColor = {0, 0, 0}, textString = "T s"), Line(points = {{-24, 0}, {54, 0}}, color = {0, 0, 0}), Line(points = {{-100, 0}, {-60, 0}}), Line(points = {{62, 0}, {100, 0}})}), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, initialScale = 0.1, grid = {2, 2}), graphics = {Line(points = {{-80, 78}, {-80, -90}}, color = {192, 192, 192}), Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 88}, {-80, 90}}), Line(points = {{-90, -80}, {82, -80}}, color = {192, 192, 192}), Polygon(lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid, points = {{90, -80}, {68, -72}, {68, -88}, {90, -80}}), Line(points = {{-80, -80}, {-80, -20}, {60, 80}}, color = {0, 0, 255}), Text(lineColor = {192, 192, 192}, extent = {{0, 6}, {60, -56}}, textString = "I"), Text(extent = {{-150, -150}, {150, -110}}, textString = "T=%T")}));
      end I;
    end Continuous;

    package Interfaces
      partial block SISO "Single Input Single Output continuous control block"
        extends CombinedCycle.Blocks.Interfaces.BlockIcon;
        CombinedCycle.Blocks.Interfaces.RealInput u "Connector of Real input signal" annotation(Placement(transformation(extent = {{-140, -20}, {-100, 20}}, rotation = 0)));
        CombinedCycle.Blocks.Interfaces.RealOutput y "Connector of Real output signal" annotation(Placement(transformation(extent = {{100, -10}, {120, 10}}, rotation = 0)));
        annotation(Window(x = 0.32, y = 0.07, width = 0.6, height = 0.6), Documentation(info = "<html>
<p>
Block has one continuous Real input and one continuous Real output signal.
</p>
</html>"));
      end SISO;

      connector RealInput = input CombinedCycle.Blocks.Interfaces.RealSignal "'input Real' as connector" annotation(defaultComponentName = "u", Icon(graphics = {Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid)}, coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = false, grid = {1, 1}, initialScale = 0.2)), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}, initialScale = 0.2), graphics = {Polygon(points = {{0, 50}, {100, 0}, {0, -50}, {0, 50}}, lineColor = {0, 0, 127}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Solid), Text(extent = {{-120, 105}, {100, 60}}, lineColor = {0, 0, 127}, textString = "%name")}), Documentation(info = "<html>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     <p>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     Connector with one input signal of type Real.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </p>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </html>"));
      connector RealOutput = output CombinedCycle.Blocks.Interfaces.RealSignal "'output Real' as connector" annotation(defaultComponentName = "y", Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics = {Polygon(points = {{-100, 100}, {100, 0}, {-100, -100}, {-100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics = {Polygon(points = {{-100, 50}, {0, 0}, {-100, -50}, {-100, 50}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-100, 140}, {130, 60}}, lineColor = {0, 0, 127}, textString = "%name")}), Documentation(info = "<html>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     <p>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     Connector with one output signal of type Real.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </p>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     </html>"));

      partial block BlockIcon "Basic graphical layout of input/output block"
        annotation(Window(x = 0, y = 0, width = 0.6, height = 0.6), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(extent = {{-100, -100}, {100, 100}}, lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-150, 150}, {150, 110}}, textString = "%name")}), Documentation(info = "<html>
<p>
Block that has only the basic icon for an input/output
block (no declarations, no equations). Most blocks
of package Modelica.Blocks inherit directly or indirectly
from this block.
</p>
</html>"));
      end BlockIcon;

      connector RealSignal "Real port (both input/output possible)"
        extends Real;
      end RealSignal;

      partial block SO "Single Output continuous control block"
        extends CombinedCycle.Blocks.Interfaces.BlockIcon;
        CombinedCycle.Blocks.Interfaces.RealOutput y "Connector of Real output signal" annotation(Placement(transformation(extent = {{100, -10}, {120, 10}}, rotation = 0)));
        annotation(Window(x = 0.25, y = 0.02, width = 0.6, height = 0.6), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics), Documentation(info = "<html>
<p>
Block has one continuous Real output signal.
</p>
</html>"));
      end SO;

      partial block SignalSource "Base class for continuous signal source"
        extends CombinedCycle.Blocks.Interfaces.SO;
        parameter Real offset = 0 "offset of output signal";
        parameter Modelica.SIunits.Time startTime = 0 "output = offset for time < startTime";
        annotation(Documentation(info = "<html>
<p>
Basic block for Real sources of package Blocks.Sources.
This component has one continuous Real output signal y
and two parameters (offset, startTime) to shift the
generated signal.
</p>
</html>"));
      end SignalSource;
    end Interfaces;

    package Sources
      block Constant "Generate constant signal of type Real"
        parameter Real k = 1 "Constant output value";
        extends CombinedCycle.Blocks.Interfaces.SO;
      equation
        y = k;
        annotation(defaultComponentName = "const", Window(x = 0.29, y = 0.19, width = 0.6, height = 0.6), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, 0}, {80, 0}}, color = {0, 0, 0}), Text(extent = {{-150, -150}, {150, -110}}, lineColor = {0, 0, 0}, textString = "k=%k")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Line(points = {{-80, 0}, {80, 0}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Text(extent = {{-75, 94}, {-22, 76}}, lineColor = {160, 160, 164}, textString = "y"), Text(extent = {{70, -80}, {94, -100}}, lineColor = {160, 160, 164}, textString = "time"), Text(extent = {{-101, 8}, {-81, -12}}, lineColor = {160, 160, 164}, textString = "k")}), Documentation(info = "<html>

</html>"));
      end Constant;

      block Step "Generate step signal of type Real"
        parameter Real height = 1 "Height of step";
        extends CombinedCycle.Blocks.Interfaces.SignalSource;
        parameter Real k = 1e5;
        constant Real pi = Modelica.Constants.pi;
      protected
        Real activation;
      equation
        activation = (atan(k * (time - startTime)) + pi / 2) / pi;
        y = offset + height * activation;
        annotation(Window(x = 0.38, y = 0.11, width = 0.6, height = 0.6), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics = {Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, -70}, {0, -70}, {0, 50}, {80, 50}}, color = {0, 0, 0}), Text(extent = {{-150, -150}, {150, -110}}, lineColor = {0, 0, 0}, textString = "startTime=%startTime")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics = {Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Line(points = {{-80, -18}, {0, -18}, {0, 50}, {80, 50}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Text(extent = {{70, -80}, {94, -100}}, lineColor = {160, 160, 164}, textString = "time"), Text(extent = {{-21, -72}, {25, -90}}, lineColor = {160, 160, 164}, textString = "startTime"), Line(points = {{0, -17}, {0, -71}}, color = {192, 192, 192}, pattern = LinePattern.Dash), Text(extent = {{-68, -36}, {-22, -54}}, lineColor = {160, 160, 164}, textString = "offset"), Line(points = {{-13, 50}, {-13, -17}}, color = {192, 192, 192}, pattern = LinePattern.Solid, thickness = 0.25, arrow = {Arrow.None, Arrow.None}), Polygon(points = {{2, 50}, {-19, 50}, {2, 50}}, lineColor = {192, 192, 192}, pattern = LinePattern.Dash), Polygon(points = {{-13, -17}, {-16, -4}, {-10, -4}, {-13, -17}, {-13, -17}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{-13, 50}, {-16, 37}, {-9, 37}, {-13, 50}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Text(extent = {{-68, 26}, {-22, 8}}, lineColor = {160, 160, 164}, textString = "height"), Polygon(points = {{-13, -69}, {-16, -56}, {-10, -56}, {-13, -69}, {-13, -69}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-13, -18}, {-13, -70}}, color = {192, 192, 192}, pattern = LinePattern.Solid, thickness = 0.25, arrow = {Arrow.None, Arrow.None}), Polygon(points = {{-13, -18}, {-16, -31}, {-9, -31}, {-13, -18}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Text(extent = {{-72, 100}, {-31, 80}}, lineColor = {160, 160, 164}, textString = "y")}), Documentation(info = "<html>

</html>"));
      end Step;

      block Ramp "Generate ramp signal"
        parameter Real height = 1 "Height of ramps";
        parameter Real duration(min = Modelica.Constants.small) = 2 "Durations of ramp";
        parameter Real offset = 0 "Offset of output signal";
        parameter Modelica.SIunits.Time startTime = 0 "Output = offset for time < startTime";
        extends CombinedCycle.Blocks.Interfaces.SO;
        parameter Real k_1 = 1e5;
        parameter Real k_2 = 1e5;
      protected
        Real activation_1;
        Real activation_2;
        constant Real pi = Modelica.Constants.pi;
      equation
        activation_1 = (atan(k_1 * (time - startTime)) + pi / 2) / pi;
        activation_2 = (atan(k_2 * (time - (startTime + duration))) + pi / 2) / pi;
        y = offset + activation_1 * (time - startTime) * height / duration * (1 - activation_2) + height * activation_2;
        annotation(Window(x = 0.19, y = 0.02, width = 0.59, height = 0.77), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics = {Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, -70}, {-40, -70}, {31, 38}}, color = {0, 0, 0}), Text(extent = {{-150, -150}, {150, -110}}, lineColor = {0, 0, 0}, textString = "duration=%duration"), Line(points = {{31, 38}, {86, 38}}, color = {0, 0, 0})}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {1, 1}), graphics = {Polygon(points = {{-80, 90}, {-88, 68}, {-72, 68}, {-80, 90}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-80, 68}, {-80, -80}}, color = {192, 192, 192}), Line(points = {{-80, -20}, {-20, -20}, {50, 50}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{-90, -70}, {82, -70}}, color = {192, 192, 192}), Polygon(points = {{90, -70}, {68, -62}, {68, -78}, {90, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{-40, -20}, {-42, -30}, {-37, -30}, {-40, -20}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Line(points = {{-40, -20}, {-40, -70}}, color = {192, 192, 192}, pattern = LinePattern.Solid, thickness = 0.25, arrow = {Arrow.None, Arrow.None}), Polygon(points = {{-40, -70}, {-43, -60}, {-38, -60}, {-40, -70}, {-40, -70}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Text(extent = {{-80, -33}, {-41, -49}}, lineColor = {160, 160, 164}, textString = "offset"), Text(extent = {{-40, -70}, {6, -88}}, lineColor = {160, 160, 164}, textString = "startTime"), Text(extent = {{-66, 92}, {-25, 72}}, lineColor = {160, 160, 164}, textString = "y"), Text(extent = {{70, -80}, {94, -100}}, lineColor = {160, 160, 164}, textString = "time"), Line(points = {{-20, -20}, {-20, -70}}, color = {192, 192, 192}, pattern = LinePattern.Dash), Line(points = {{-19, -20}, {50, -20}}, color = {192, 192, 192}, pattern = LinePattern.Solid, thickness = 0.25, arrow = {Arrow.None, Arrow.None}), Line(points = {{50, 50}, {101, 50}}, color = {0, 0, 0}, thickness = 0.5), Line(points = {{50, 50}, {50, -20}}, color = {192, 192, 192}, pattern = LinePattern.Solid, thickness = 0.25, arrow = {Arrow.None, Arrow.None}), Polygon(points = {{50, -20}, {42, -18}, {42, -22}, {50, -20}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{-20, -20}, {-11, -18}, {-11, -22}, {-20, -20}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{50, 50}, {48, 40}, {53, 40}, {50, 50}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{50, -20}, {47, -10}, {52, -10}, {50, -20}, {50, -20}}, lineColor = {192, 192, 192}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Text(extent = {{53, 25}, {82, 7}}, lineColor = {160, 160, 164}, textString = "height"), Text(extent = {{0, -17}, {35, -37}}, lineColor = {160, 160, 164}, textString = "duration")}), Documentation(info = "<html>

</html>"));
      end Ramp;

      model BandLimitedWhiteNoise "Band-limited, zero-mean white noise generator (algorithm-based)"
        extends Modelica.Blocks.Interfaces.SO;
        parameter Modelica.SIunits.AngularFrequency omega_c "Noise bandwidth";
        parameter Real A "Mean amplitude of noise signal";
        final parameter Integer N = 18 "Number of bits in PRBS generator";
        parameter Boolean seed[N] = {false, false, false, true, true, false, true, false, true, true, true, false, true, false, false, true, false, true} "Seed of PRBS generator";
        final parameter Modelica.SIunits.Time Tc = 0.3 / omega_c "Sampling time for PRBS generator";
        final parameter Real M = 1 / sqrt(0.3) "Scaling factor";
        discrete Real yPRBS "Segnale PRBS";
        Real x[2](each start = 0, each fixed = true) "Low-pass filter state";
        Boolean s[N] "Shift register";
        Boolean h "Temporary variable";
      initial algorithm
        for i in 1:N loop
          s[i] := seed[i];
        end for;
      algorithm
        when sample(-50, Tc) then
          yPRBS := if s[1] then M * A else -M * A;
          h := not s[18] and s[5] or s[18] and not s[5];
          h := not s[2] and h or s[2] and not h;
          h := not s[1] and h or s[1] and not h;
          for j in N:(-1):2 loop
            s[j] := s[j - 1];
          end for;
          s[1] := h;
        end when;
      equation
        der(x[1]) = x[2];
        der(x[2]) / omega_c ^ 2 + sqrt(2) / omega_c * x[2] + x[1] = yPRBS;
        y = x[1];
      end BandLimitedWhiteNoise;

      model BandLimitedWhiteNoiseEquations "Band-limited, zero-mean white noise generator (equation-based)"
        extends Modelica.Blocks.Interfaces.SO;
        parameter Modelica.SIunits.AngularFrequency omega_c "Noise bandwidth";
        parameter Real A "Mean amplitude of noise signal";
        final parameter Integer N = 18 "Number of bits in PRBS generator";
        parameter Boolean seed[N] = {false, false, false, true, true, false, true, false, true, true, true, false, true, false, false, true, false, true} "Seed of PRBS generator";
        final parameter Modelica.SIunits.Time Tc = 0.3 / omega_c "Sampling time for PRBS generator";
        final parameter Real M = 1 / sqrt(0.3) "Scaling factor";
        discrete Real yPRBS "Segnale PRBS";
        Real x[2](each start = 0, each fixed = true) "Low-pass filter state";
        Boolean s[N](start = seed, fixed = true) "Shift register";
        Boolean h[3] "Temporary variable";
      equation
        when sample(-50, Tc) then
          yPRBS = if pre(s[1]) then M * A else -M * A;
          h[1] = not pre(s[18]) and pre(s[5]) or pre(s[18]) and not pre(s[5]);
          h[2] = not pre(s[2]) and h[1] or pre(s[2]) and not h[1];
          h[3] = not pre(s[1]) and h[2] or pre(s[1]) and not h[2];
          for j in N:(-1):2 loop
            s[j] = pre(s[j - 1]);
          end for;
          s[1] = h[3];
        end when;
        der(x[1]) = x[2];
        der(x[2]) / omega_c ^ 2 + sqrt(2) / omega_c * x[2] + x[1] = yPRBS;
        y = x[1];
      end BandLimitedWhiteNoiseEquations;

      package Test
        model TestNoiseGenerators
          BandLimitedWhiteNoise noise1(omega_c = 0.2, A = 0.01) annotation(Placement(transformation(extent = {{-20, 20}, {0, 40}})));
          BandLimitedWhiteNoiseEquations noise2(omega_c = 0.2, A = 0.01) annotation(Placement(transformation(extent = {{-20, -20}, {0, 0}})));
          annotation(experiment(StopTime = 1000, Tolerance = 1e-006), __Dymola_experimentSetupOutput);
        end TestNoiseGenerators;
      end Test;
    end Sources;

    package Math
      block Feedback "Output difference between commanded and feedback input"
        input Modelica.Blocks.Interfaces.RealInput u1 annotation(Placement(transformation(extent = {{-100, -20}, {-60, 20}}, rotation = 0)));
        input Modelica.Blocks.Interfaces.RealInput u2 annotation(Placement(transformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
        output Modelica.Blocks.Interfaces.RealOutput y annotation(Placement(transformation(extent = {{80, -10}, {100, 10}}, rotation = 0)));
      equation
        y = u1 - u2;
        annotation(Window(x = 0.35, y = 0.02, width = 0.52, height = 0.68), Documentation(info = "
<HTML>
<p>
This blocks computes output <b>y</b> as <i>difference</i> of the
commanded input <b>u1</b> and the feedback
input <b>u2</b>:
</p>
<pre>
    <b>y</b> = <b>u1</b> - <b>u2</b>;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   n = 2
  results in the following equations:
     y = u1 - u2
</pre>

</HTML>
"), Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Ellipse(extent = {{-20, 20}, {20, -20}}, lineColor = {0, 0, 127}, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid), Line(points = {{-60, 0}, {-20, 0}}, color = {0, 0, 127}), Line(points = {{20, 0}, {80, 0}}, color = {0, 0, 127}), Line(points = {{0, -20}, {0, -60}}, color = {0, 0, 127}), Text(extent = {{-14, 0}, {82, -94}}, lineColor = {0, 0, 0}, textString = "-"), Text(extent = {{-100, 110}, {100, 60}}, textString = "%name")}), Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}, grid = {2, 2}), graphics = {Ellipse(extent = {{-20, 20}, {20, -20}}, lineColor = {0, 0, 255}, pattern = LinePattern.Solid, lineThickness = 0.25, fillColor = {235, 235, 235}, fillPattern = FillPattern.Solid), Line(points = {{-60, 0}, {-20, 0}}), Line(points = {{20, 0}, {80, 0}}), Line(points = {{0, -20}, {0, -60}}), Text(extent = {{-12, 10}, {84, -84}}, lineColor = {0, 0, 0}, textString = "-")}));
      end Feedback;
    end Math;

    model Pseudo_controller
      CombinedCycle.Blocks.Interfaces.RealInput T0 "Desired Value for Tit" annotation(Placement(visible = true, transformation(origin = {-78, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 1.11022e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      CombinedCycle.Blocks.Interfaces.RealInput Tm "Measured Value for Tit" annotation(Placement(visible = true, transformation(origin = {-6, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 80}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      CombinedCycle.Blocks.Interfaces.RealOutput w(min = 0) "Resulting Flow Rate" annotation(Placement(visible = true, transformation(origin = {70, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      parameter Boolean Attemperation = true "Select false to deactivate the attemperation";
    equation
      if Attemperation == true then
        Tm = T0;
      else
        w = 0;
      end if;
      annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(origin = {-20, 20}, lineColor = {0, 170, 0}, fillColor = {0, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-40, 40}, {80, -80}})}));
    end Pseudo_controller;
  end Blocks;

  package Circuits_Tests
    model CombinedCycle_ThermalStress
      Components.Gas.BraytonJoule_exhaustGases braytonJoule_exhaustGases(T_min = 548, T_max = 843, w_min = 454, w_max = 585) annotation(Placement(transformation(extent = {{-60, 62}, {-40, 82}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst_1(Tstart = 600 + 273.15, Tmstart = 600 + 273.15, V = 10) annotation(Placement(transformation(extent = {{-30, 82}, {-10, 62}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst_2(V = 10, Tstart = 500 + 273.15, Tmstart = 500 + 273.15) annotation(Placement(transformation(extent = {{6, 82}, {26, 62}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst_3(V = 10, Tstart = 400 + 273.15, Tmstart = 400 + 273.15) annotation(Placement(transformation(extent = {{40, 82}, {60, 62}}, rotation = 0)));
      Components.Gas.SinkP environment annotation(Placement(transformation(extent = {{76, 62}, {96, 82}}, rotation = 0)));
      Components.Water.Flow1DLiquidPConstWall Economizer(pstart = 90e5, Tstart = 400, Tmstart = 400, wnom = 70, V = 30, Cm = 200 * Economizer.wnom / 4186 / 50) annotation(Placement(transformation(extent = {{58, -52}, {38, -32}}, rotation = 0)));
      Components.Water.Drum Evaporator(V = 15, pstart = 60e5, alphastart = 0.5, Cm = 200 * 70 / 4186 / 5, Tstart = 600, Tmstart = 600) annotation(Placement(transformation(extent = {{26, -50}, {6, -30}}, rotation = 0)));
      Components.Water.TurbineThermalPort turbine(pstart = 70e5, h_iso = 20e5, Gamma = 100, k = 70 / 90e5) annotation(Placement(transformation(extent = {{-70, -52}, {-50, -32}}, rotation = 0)));
      Components.Thermal.ThermalConductor thermalConductor_3(G = 350000) annotation(Placement(transformation(origin = {50, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.ThermalConductor thermalConductor_1(G = 500000) annotation(Placement(transformation(origin = {-20, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Blocks.Sources.Ramp load(offset = 0.5, height = -0.3, startTime = 1000, duration = 2000) annotation(Placement(transformation(extent = {{-92, 62}, {-72, 82}}, rotation = 0)));
      Components.Water.Flow1DVapourPConstWall SuperHeater(pstart = 90e5, Tstart = 800, Tmstart = 800, V = 5, Cm = 200 * 70 / 4186 / 5) annotation(Placement(transformation(extent = {{-10, -52}, {-30, -32}}, rotation = 0)));
      Components.Thermal.ThermalConductor thermalConductor_2(G = 1300000) annotation(Placement(transformation(origin = {16, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.CylinderThermalStress_radial_noVector cylinderThermalStress(rint = 0.4, rext = 0.6, Tstartint = 800, Tstartext = 800) annotation(Placement(transformation(extent = {{-70, -70}, {-50, -50}}, rotation = 0)));
      Blocks.Continuous.PI PI(y_start = 70, k = -1000, T = 0.01) annotation(Placement(transformation(extent = {{-42, -6}, {-24, 10}}, rotation = 0)));
      Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{-68, -8}, {-48, 12}}, rotation = 0)));
      Blocks.Sources.Step alpha_SP(height = 0, startTime = 10, offset = 0.5) annotation(Placement(transformation(extent = {{-94, -8}, {-74, 12}}, rotation = 0)));
      Components.Water.SourceW_waterLiquid_Input WaterLiquidSource(wstart = 70, Tnom = 273.15 + 35) annotation(Placement(transformation(extent = {{90, -52}, {70, -32}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow_DHT_N1 fixedHeatFlow_DHT_N1_1(Q_flow = 0) annotation(Placement(transformation(origin = {-60, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(flow1DGasPConst_3.outlet, environment.flange) annotation(Line(points = {{60, 72}, {76, 72}}, color = {159, 159, 223}));
      connect(flow1DGasPConst_2.outlet, flow1DGasPConst_3.inlet) annotation(Line(points = {{26, 72}, {40, 72}}, color = {159, 159, 223}));
      connect(flow1DGasPConst_1.outlet, flow1DGasPConst_2.inlet) annotation(Line(points = {{-10, 72}, {6, 72}}, color = {159, 159, 223}));
      connect(braytonJoule_exhaustGases.flange, flow1DGasPConst_1.inlet) annotation(Line(points = {{-40, 72}, {-30, 72}}, color = {159, 159, 223}));
      connect(thermalConductor_3.port_b, flow1DGasPConst_3.wall) annotation(Line(points = {{50, 52}, {50, 67}}, color = {191, 0, 0}));
      connect(thermalConductor_3.port_a, Economizer.wall) annotation(Line(points = {{50, 32}, {50, -37}, {48, -37}}, color = {191, 0, 0}));
      connect(thermalConductor_1.port_b, flow1DGasPConst_1.wall) annotation(Line(points = {{-20, 52}, {-20, 67}}, color = {191, 0, 0}));
      connect(flow1DGasPConst_2.wall, thermalConductor_2.port_b) annotation(Line(points = {{16, 67}, {16, 52}}, color = {191, 0, 0}));
      connect(cylinderThermalStress.internalBoundary, turbine.dHT) annotation(Line(points = {{-60, -57}, {-60, -46.2}}, color = {255, 127, 0}));
      connect(load.y, braytonJoule_exhaustGases.load) annotation(Line(points = {{-71, 72}, {-60, 72}}, color = {0, 0, 127}));
      connect(Evaporator.inlet, Economizer.outlet) annotation(Line(points = {{24, -45.8}, {30.7, -45.8}, {30.7, -42}, {38, -42}}, color = {159, 159, 223}));
      connect(Evaporator.outlet, SuperHeater.inlet) annotation(Line(points = {{10, -34}, {0.4, -34}, {0.4, -42}, {-10, -42}}, color = {159, 159, 223}));
      connect(SuperHeater.wall, thermalConductor_1.port_a) annotation(Line(points = {{-20, -37}, {-20, 32}}, color = {191, 0, 0}));
      connect(Evaporator.wall, thermalConductor_2.port_a) annotation(Line(points = {{16, -48.9}, {16, 32}}, color = {191, 0, 0}));
      connect(turbine.inlet, SuperHeater.outlet) annotation(Line(points = {{-51, -42}, {-30, -42}}, color = {159, 159, 223}));
      connect(feedback.y, PI.u) annotation(Line(points = {{-49, 2}, {-46.4, 2}, {-43.8, 2}}, color = {0, 0, 127}));
      connect(alpha_SP.y, feedback.u1) annotation(Line(points = {{-73, 2}, {-71.25, 2}, {-71.25, 2}, {-69.5, 2}, {-69.5, 2}, {-66, 2}}, color = {0, 0, 127}));
      connect(Evaporator.voidFraction, feedback.u2) annotation(Line(points = {{8.4, -46}, {-6, -46}, {-6, -20}, {-58, -20}, {-58, -6}}, color = {0, 0, 127}));
      connect(PI.y, WaterLiquidSource.w) annotation(Line(points = {{-23.1, 2}, {86.2, 2}, {86.2, -36}}, color = {0, 0, 127}));
      connect(Economizer.inlet, WaterLiquidSource.flange) annotation(Line(points = {{58, -42}, {70, -42}}, color = {159, 159, 223}));
      connect(fixedHeatFlow_DHT_N1_1.port, cylinderThermalStress.externalBoundary) annotation(Line(points = {{-60, -72}, {-60, -63}}, color = {255, 127, 0}));
      annotation(Diagram(graphics), experiment(StopTime = 1000), experimentSetupOutput);
    end CombinedCycle_ThermalStress;

    model CombinedCycle_ThermalStress_startUp
      Components.Gas.BraytonJoule_exhaustGases braytonJoule_exhaustGases(T_min = 548, T_max = 843, w_min = 454, w_max = 585, wstart = 460, Tstart = 550) annotation(Placement(transformation(extent = {{-60, 68}, {-40, 88}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst_1(V = 10, Tstart = 625, Tmstart = 625) annotation(Placement(transformation(extent = {{-30, 88}, {-10, 68}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst_2(V = 10, Tstart = 500, Tmstart = 500) annotation(Placement(transformation(extent = {{6, 88}, {26, 68}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst_3(V = 10, Tstart = 450, Tmstart = 450) annotation(Placement(transformation(extent = {{40, 88}, {60, 68}}, rotation = 0)));
      Components.Gas.SinkP environment annotation(Placement(transformation(extent = {{76, 68}, {96, 88}}, rotation = 0)));
      Components.Water.Flow1DLiquidPConstWall Economizer(V = 30, Cm = 200 * 70 / 4186 / 50, pstart = 3e6, Tstart = 400, Tmstart = 400, wnom = 30) annotation(Placement(transformation(extent = {{58, -48}, {38, -28}}, rotation = 0)));
      Components.Water.Drum Evaporator(V = 15, alphastart = 0.5, Cm = 200 * 70 / 4186 / 5, pstart = 3e6, Tstart = 450, Tmstart = 450) annotation(Placement(transformation(extent = {{26, -46}, {6, -26}}, rotation = 0)));
      Components.Water.TurbineThermalPort turbine(h_iso = 20e5, Gamma = 100, k = 70 / 90e5, pstart = 3e6) annotation(Placement(transformation(extent = {{-70, -46}, {-50, -26}}, rotation = 0)));
      Components.Thermal.ThermalConductor thermalConductor_3(G = 270000) annotation(Placement(transformation(origin = {50, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.ThermalConductor thermalConductor_1(G = 500000) annotation(Placement(transformation(origin = {-20, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Blocks.Sources.Ramp load(offset = 0.11, height = 0.88, startTime = 1, duration = 10000) annotation(Placement(transformation(extent = {{-92, 68}, {-72, 88}}, rotation = 0)));
      Components.Water.Flow1DVapourPConstWall SuperHeater(V = 5, Cm = 200 * 70 / 4186 / 5, pstart = 3e6, Tstart = 550, Tmstart = 550) annotation(Placement(transformation(extent = {{-10, -46}, {-30, -26}}, rotation = 0)));
      Components.Thermal.ThermalConductor thermalConductor_2(G = 1300000) annotation(Placement(transformation(origin = {16, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.CylinderThermalStress_radial_noVector cylinderThermalStress(Tstartint = 550, Tstartext = 550, rint = 0.04, rext = 0.3) annotation(Placement(transformation(extent = {{-70, -68}, {-50, -48}}, rotation = 0)));
      Blocks.Continuous.PI PI(k = -1000, T = 0.01, y_start = 30) annotation(Placement(transformation(extent = {{-42, 0}, {-24, 16}}, rotation = 0)));
      Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{-68, -2}, {-48, 18}}, rotation = 0)));
      Blocks.Sources.Step alpha_SP(height = 0, offset = 0.5, startTime = 0) annotation(Placement(transformation(extent = {{-94, -2}, {-74, 18}}, rotation = 0)));
      Components.Water.SourceW_waterLiquid_Input WaterLiquidSource(Tnom = 273.15 + 35, pstart = 3e6, wstart = 30) annotation(Placement(transformation(extent = {{90, -48}, {70, -28}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow_DHT_N1 fixedHeatFlow_DHT_N1_1(Q_flow = 0) annotation(Placement(transformation(origin = {-60, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    initial equation
      der(cylinderThermalStress.T__2) = 0;
      der(Economizer.wat_liq_out.T) = 0;
      der(flow1DGasPConst_1.gas_out.T) = 0;
      der(flow1DGasPConst_2.gas_out.T) = 0;
      der(flow1DGasPConst_3.gas_out.T) = 0;
      der(SuperHeater.wat_vap_out.T) = 0;
    equation
      connect(flow1DGasPConst_3.outlet, environment.flange) annotation(Line(points = {{60, 78}, {76, 78}}, color = {159, 159, 223}));
      connect(flow1DGasPConst_2.outlet, flow1DGasPConst_3.inlet) annotation(Line(points = {{26, 78}, {40, 78}}, color = {159, 159, 223}));
      connect(flow1DGasPConst_1.outlet, flow1DGasPConst_2.inlet) annotation(Line(points = {{-10, 78}, {6, 78}}, color = {159, 159, 223}));
      connect(braytonJoule_exhaustGases.flange, flow1DGasPConst_1.inlet) annotation(Line(points = {{-40, 78}, {-30, 78}}, color = {159, 159, 223}));
      connect(thermalConductor_3.port_b, flow1DGasPConst_3.wall) annotation(Line(points = {{50, 58}, {50, 73}}, color = {191, 0, 0}));
      connect(thermalConductor_3.port_a, Economizer.wall) annotation(Line(points = {{50, 38}, {50, -33}, {48, -33}}, color = {191, 0, 0}));
      connect(thermalConductor_1.port_b, flow1DGasPConst_1.wall) annotation(Line(points = {{-20, 58}, {-20, 73}}, color = {191, 0, 0}));
      connect(flow1DGasPConst_2.wall, thermalConductor_2.port_b) annotation(Line(points = {{16, 73}, {16, 58}}, color = {191, 0, 0}));
      connect(cylinderThermalStress.internalBoundary, turbine.dHT) annotation(Line(points = {{-60, -55}, {-60, -40.2}}, color = {255, 127, 0}));
      connect(load.y, braytonJoule_exhaustGases.load) annotation(Line(points = {{-71, 78}, {-60, 78}}, color = {0, 0, 127}));
      connect(Evaporator.inlet, Economizer.outlet) annotation(Line(points = {{24, -41.8}, {30.7, -41.8}, {30.7, -38}, {38, -38}}, color = {159, 159, 223}));
      connect(Evaporator.outlet, SuperHeater.inlet) annotation(Line(points = {{10, -30}, {0.4, -30}, {0.4, -36}, {-10, -36}}, color = {159, 159, 223}));
      connect(SuperHeater.wall, thermalConductor_1.port_a) annotation(Line(points = {{-20, -31}, {-20, 38}}, color = {191, 0, 0}));
      connect(Evaporator.wall, thermalConductor_2.port_a) annotation(Line(points = {{16, -44.9}, {16, 38}}, color = {191, 0, 0}));
      connect(turbine.inlet, SuperHeater.outlet) annotation(Line(points = {{-51, -36}, {-30, -36}}, color = {159, 159, 223}));
      connect(feedback.y, PI.u) annotation(Line(points = {{-49, 8}, {-43.8, 8}}, color = {0, 0, 127}));
      connect(alpha_SP.y, feedback.u1) annotation(Line(points = {{-73, 8}, {-66, 8}}, color = {0, 0, 127}));
      connect(Evaporator.voidFraction, feedback.u2) annotation(Line(points = {{8.4, -42}, {-6, -42}, {-6, -14}, {-58, -14}, {-58, 0}}, color = {0, 0, 127}));
      connect(PI.y, WaterLiquidSource.w) annotation(Line(points = {{-23.1, 8}, {86.2, 8}, {86.2, -32}}, color = {0, 0, 127}));
      connect(Economizer.inlet, WaterLiquidSource.flange) annotation(Line(points = {{58, -38}, {70, -38}}, color = {159, 159, 223}));
      connect(fixedHeatFlow_DHT_N1_1.port, cylinderThermalStress.externalBoundary) annotation(Line(points = {{-60, -76}, {-60, -61}}, color = {255, 127, 0}));
      annotation(Diagram(graphics), experiment(StopTime = 1000), experimentSetupOutput);
    end CombinedCycle_ThermalStress_startUp;

    model CombinedCycle_sub1
      Components.Gas.BraytonJoule_exhaustGases braytonJoule_exhaustGases(T_min = 600, T_max = 800, w_min = 300, w_max = 500) annotation(Placement(transformation(extent = {{-64, 24}, {-44, 44}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst(V = 5) annotation(Placement(transformation(extent = {{-36, 44}, {-16, 24}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst1(V = 5) annotation(Placement(transformation(extent = {{-2, 44}, {18, 24}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst2(V = 5) annotation(Placement(transformation(extent = {{34, 44}, {54, 24}}, rotation = 0)));
      Components.Gas.SinkP sinkP annotation(Placement(transformation(extent = {{74, 24}, {94, 44}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant ramp(k = 0.5) annotation(Placement(transformation(extent = {{-96, 24}, {-76, 44}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 0) annotation(Placement(transformation(origin = {-26, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow1(Q_flow = 0) annotation(Placement(transformation(origin = {8, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow2(Q_flow = 0) annotation(Placement(transformation(origin = {44, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(flow1DGasPConst2.outlet, sinkP.flange) annotation(Line(points = {{54, 34}, {74, 34}}, color = {159, 159, 223}));
      connect(flow1DGasPConst1.outlet, flow1DGasPConst2.inlet) annotation(Line(points = {{18, 34}, {34, 34}}, color = {159, 159, 223}));
      connect(flow1DGasPConst.outlet, flow1DGasPConst1.inlet) annotation(Line(points = {{-16, 34}, {-2, 34}}, color = {159, 159, 223}));
      connect(braytonJoule_exhaustGases.flange, flow1DGasPConst.inlet) annotation(Line(points = {{-44, 34}, {-36, 34}}, color = {159, 159, 223}));
      connect(ramp.y, braytonJoule_exhaustGases.load) annotation(Line(points = {{-75, 34}, {-64, 34}}, color = {0, 0, 127}));
      connect(fixedHeatFlow2.port, flow1DGasPConst2.wall) annotation(Line(points = {{44, 6}, {44, 29}}, color = {191, 0, 0}));
      connect(flow1DGasPConst.wall, fixedHeatFlow.port) annotation(Line(points = {{-26, 29}, {-26, 8}, {-26, 8}}, color = {191, 0, 0}));
      connect(fixedHeatFlow1.port, flow1DGasPConst1.wall) annotation(Line(points = {{8, 6}, {8, 17.5}, {8, 17.5}, {8, 29}}, color = {191, 0, 0}));
      annotation(Diagram(graphics));
    end CombinedCycle_sub1;

    model CombinedCycle_sub2
      Components.Water.Flow1DLiquidPConstWall Economizer(pstart = 90e5, Tstart = 400, Tmstart = 400, wnom = 70, V = 30, Cm = 200 * Economizer.wnom / 4186 / 50) annotation(Placement(transformation(extent = {{50, -50}, {30, -30}}, rotation = 0)));
      Components.Water.Flow1DVapourPConstWall SuperHeater(pstart = 90e5, Tstart = 800, Tmstart = 800, V = 5, Cm = 200 * 70 / 4186 / 5) annotation(Placement(transformation(extent = {{-32, -50}, {-52, -30}}, rotation = 0)));
      Components.Water.Turbine turbine(k = 10, h_iso = 2.3e5) annotation(Placement(transformation(extent = {{-86, -50}, {-66, -30}}, rotation = 0)));
      Components.Water.SourceW_waterLiquid_Input sourceW(Tnom = 273.15 + 35, wstart = 70) annotation(Placement(transformation(extent = {{88, -50}, {68, -30}}, rotation = 0)));
      Components.Water.Drum Evaporator(V = 15, pstart = 60e5, alphastart = 0.5, Cm = 200 * 70 / 4186 / 5, Tstart = 600, Tmstart = 600) annotation(Placement(transformation(extent = {{6, -50}, {-14, -30}}, rotation = 0)));
      Blocks.Continuous.PI PI(k = -1000, T = 0.01, y_start = 30) annotation(Placement(transformation(extent = {{-44, 30}, {-26, 46}}, rotation = 0)));
      Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{-70, 28}, {-50, 48}}, rotation = 0)));
      Blocks.Sources.Step alpha_SP(height = 0, offset = 0.5, startTime = 0) annotation(Placement(transformation(extent = {{-96, 28}, {-76, 48}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 3e7) annotation(Placement(transformation(origin = {40, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow1(Q_flow = 9e7) annotation(Placement(transformation(origin = {-4, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow2(Q_flow = 4e7) annotation(Placement(transformation(origin = {-42, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    equation
      connect(turbine.inlet, SuperHeater.outlet) annotation(Line(points = {{-67, -40}, {-52, -40}}, color = {159, 159, 223}));
      connect(Economizer.inlet, sourceW.flange) annotation(Line(points = {{50, -40}, {68, -40}}, color = {159, 159, 223}));
      connect(SuperHeater.inlet, Evaporator.outlet) annotation(Line(points = {{-32, -40}, {-22, -40}, {-22, -34.2}, {-11.2, -34.2}}, color = {159, 159, 223}));
      connect(Evaporator.inlet, Economizer.outlet) annotation(Line(points = {{3.4, -45.4}, {15.7, -45.4}, {15.7, -40}, {30, -40}}, color = {159, 159, 223}));
      connect(Evaporator.voidFraction, feedback.u2) annotation(Line(points = {{-11.6, -46}, {-24, -46}, {-24, -16}, {-60, -16}, {-60, 30}}, color = {0, 0, 127}));
      connect(PI.y, sourceW.w) annotation(Line(points = {{-25.1, 38}, {84.2, 38}, {84.2, -34}}, color = {0, 0, 127}));
      connect(fixedHeatFlow.port, Economizer.wall) annotation(Line(points = {{40, -24}, {40, -35}}, color = {191, 0, 0}));
      connect(Evaporator.wall, fixedHeatFlow1.port) annotation(Line(points = {{-4, -31}, {-4, -28.75}, {-4, -28.75}, {-4, -26.5}, {-4, -22}, {-4, -22}}, color = {191, 0, 0}));
      connect(alpha_SP.y, feedback.u1) annotation(Line(points = {{-75, 38}, {-68, 38}}, color = {0, 0, 127}));
      connect(feedback.y, PI.u) annotation(Line(points = {{-51, 38}, {-45.8, 38}}, color = {0, 0, 127}));
      connect(SuperHeater.wall, fixedHeatFlow2.port) annotation(Line(points = {{-42, -35}, {-42, -20}}, color = {191, 0, 0}));
      annotation(Diagram(graphics));
    end CombinedCycle_sub2;

    model CombinedCycle_sub3
      Components.Water.Flow1DLiquidPConstWall Economizer(pstart = 90e5, Tstart = 400, Tmstart = 400, wnom = 70, V = 30, Cm = 200 * Economizer.wnom / 4186 / 50) annotation(Placement(transformation(extent = {{50, -50}, {30, -30}}, rotation = 0)));
      Components.Water.Flow1DVapourPConstWall SuperHeater(pstart = 90e5, Tstart = 800, Tmstart = 800, V = 5, Cm = 200 * 70 / 4186 / 5) annotation(Placement(transformation(extent = {{-32, -50}, {-52, -30}}, rotation = 0)));
      Components.Water.SourceW_waterLiquid_Input sourceW(Tnom = 273.15 + 35, wstart = 70) annotation(Placement(transformation(extent = {{88, -50}, {68, -30}}, rotation = 0)));
      Components.Water.Drum Evaporator(V = 15, pstart = 60e5, alphastart = 0.5, Cm = 200 * 70 / 4186 / 5, Tstart = 600, Tmstart = 600) annotation(Placement(transformation(extent = {{6, -50}, {-14, -30}}, rotation = 0)));
      Blocks.Continuous.PI PI(k = -1000, T = 0.01, y_start = 30) annotation(Placement(transformation(extent = {{-44, 30}, {-26, 46}}, rotation = 0)));
      Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{-70, 28}, {-50, 48}}, rotation = 0)));
      Blocks.Sources.Step alpha_SP(height = 0, offset = 0.5, startTime = 0) annotation(Placement(transformation(extent = {{-96, 28}, {-76, 48}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 3e7) annotation(Placement(transformation(origin = {40, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow1(Q_flow = 9e7) annotation(Placement(transformation(origin = {-4, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow2(Q_flow = 4e7) annotation(Placement(transformation(origin = {-42, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Water.TurbineThermalPort turbine1(h_iso = 20e5, Gamma = 100, k = 70 / 90e5, pstart = 3e6) annotation(Placement(transformation(extent = {{-86, -50}, {-66, -30}}, rotation = 0)));
      Components.Thermal.CylinderThermalStress_radial_noVector cylinderThermalStress(Tstartint = 550, Tstartext = 550, rint = 0.04, rext = 0.3) annotation(Placement(transformation(extent = {{-86, -74}, {-66, -54}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow_DHT_N1 fixedHeatFlow_DHT_N1_1(Q_flow = 0) annotation(Placement(transformation(origin = {-76, -88}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(Economizer.inlet, sourceW.flange) annotation(Line(points = {{50, -40}, {68, -40}}, color = {159, 159, 223}));
      connect(SuperHeater.inlet, Evaporator.outlet) annotation(Line(points = {{-32, -40}, {-22, -40}, {-22, -34.2}, {-11.2, -34.2}}, color = {159, 159, 223}));
      connect(Evaporator.inlet, Economizer.outlet) annotation(Line(points = {{3.4, -45.4}, {15.7, -45.4}, {15.7, -40}, {30, -40}}, color = {159, 159, 223}));
      connect(Evaporator.voidFraction, feedback.u2) annotation(Line(points = {{-11.6, -46}, {-24, -46}, {-24, -16}, {-60, -16}, {-60, 30}}, color = {0, 0, 127}));
      connect(PI.y, sourceW.w) annotation(Line(points = {{-25.1, 38}, {84.2, 38}, {84.2, -34}}, color = {0, 0, 127}));
      connect(fixedHeatFlow.port, Economizer.wall) annotation(Line(points = {{40, -24}, {40, -35}}, color = {191, 0, 0}));
      connect(Evaporator.wall, fixedHeatFlow1.port) annotation(Line(points = {{-4, -31}, {-4, -28.75}, {-4, -28.75}, {-4, -26.5}, {-4, -22}, {-4, -22}}, color = {191, 0, 0}));
      connect(alpha_SP.y, feedback.u1) annotation(Line(points = {{-75, 38}, {-68, 38}}, color = {0, 0, 127}));
      connect(feedback.y, PI.u) annotation(Line(points = {{-51, 38}, {-45.8, 38}}, color = {0, 0, 127}));
      connect(SuperHeater.wall, fixedHeatFlow2.port) annotation(Line(points = {{-42, -35}, {-42, -20}}, color = {191, 0, 0}));
      connect(cylinderThermalStress.internalBoundary, turbine1.dHT) annotation(Line(points = {{-76, -61}, {-76, -44.2}}, color = {255, 127, 0}));
      connect(turbine1.inlet, SuperHeater.outlet) annotation(Line(points = {{-67, -40}, {-52, -40}}, color = {159, 159, 223}));
      connect(fixedHeatFlow_DHT_N1_1.port, cylinderThermalStress.externalBoundary) annotation(Line(points = {{-76, -78}, {-76, -67}}, color = {255, 127, 0}));
      annotation(Diagram(graphics));
    end CombinedCycle_sub3;

    model GasFlow1d_test
      Components.Gas.SourceP sourceP(pnom = 120000) annotation(Placement(transformation(extent = {{-80, 16}, {-60, 36}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst(V = 1) annotation(Placement(transformation(extent = {{-34, 36}, {-14, 16}}, rotation = 0)));
      Components.Gas.SinkP sinkP annotation(Placement(transformation(extent = {{80, 16}, {100, 36}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst1(V = 1) annotation(Placement(transformation(extent = {{6, 36}, {26, 16}}, rotation = 0)));
      Components.Gas.Flow1DGasPConst flow1DGasPConst2(V = 1) annotation(Placement(transformation(extent = {{42, 36}, {62, 16}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 0) annotation(Placement(transformation(origin = {-24, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow1(Q_flow = 0) annotation(Placement(transformation(origin = {16, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow2(Q_flow = 0) annotation(Placement(transformation(origin = {52, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    equation
      connect(sourceP.flange, flow1DGasPConst.inlet) annotation(Line(points = {{-60, 26}, {-34, 26}}, color = {159, 159, 223}));
      connect(flow1DGasPConst.outlet, flow1DGasPConst1.inlet) annotation(Line(points = {{-14, 26}, {6, 26}}, color = {159, 159, 223}));
      connect(flow1DGasPConst1.outlet, flow1DGasPConst2.inlet) annotation(Line(points = {{26, 26}, {42, 26}}, color = {159, 159, 223}));
      connect(flow1DGasPConst2.outlet, sinkP.flange) annotation(Line(points = {{62, 26}, {80, 26}}, color = {159, 159, 223}));
      connect(fixedHeatFlow.port, flow1DGasPConst.wall) annotation(Line(points = {{-24, 5.55112e-16}, {-24, 21}}, color = {191, 0, 0}));
      connect(fixedHeatFlow1.port, flow1DGasPConst1.wall) annotation(Line(points = {{16, 5.55112e-16}, {16, 21}}, color = {191, 0, 0}));
      connect(fixedHeatFlow2.port, flow1DGasPConst2.wall) annotation(Line(points = {{52, 5.55112e-16}, {52, 21}}, color = {191, 0, 0}));
      annotation(Diagram(graphics));
    end GasFlow1d_test;

    model WaterLiquidFlow1d_test
      Components.Gas.SinkP sinkP annotation(Placement(transformation(origin = {-56, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Components.Water.Flow1DLiquidPConstWall Economizer(V = 1) annotation(Placement(transformation(extent = {{6, -20}, {-14, 0}}, rotation = 0)));
      Components.Water.SourceW_waterLiquid sourceW(wnom = 1) annotation(Placement(transformation(extent = {{56, -20}, {36, 0}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 0) annotation(Placement(transformation(origin = {-4, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    equation
      connect(Economizer.inlet, sourceW.flange) annotation(Line(points = {{6, -10}, {36, -10}}, color = {159, 159, 223}));
      connect(sinkP.flange, Economizer.outlet) annotation(Line(points = {{-46, -10}, {-30, -10}, {-30, -10}, {-14, -10}}, color = {159, 159, 223}));
      connect(fixedHeatFlow.port, Economizer.wall) annotation(Line(points = {{-4, 12}, {-4, 7.75}, {-4, 7.75}, {-4, 3.5}, {-4, -5}, {-4, -5}}, color = {191, 0, 0}));
      annotation(Diagram(graphics));
    end WaterLiquidFlow1d_test;

    model WaterDrum_test
      Components.Water.Drum Evaporator(pstart = 70e5, Tstart = 800, Tmstart = 800, alphastart = 0.3, V = 30, Cm = 2000 * 4186) annotation(Placement(transformation(extent = {{18, -38}, {-2, -18}}, rotation = 0)));
      Components.Water.SourceW_waterVapour_Input sourceW(wstart = 70, Tnom = 700) annotation(Placement(transformation(extent = {{70, -40}, {50, -20}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 1.2e7) annotation(Placement(transformation(origin = {8, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
      Components.Water.SinkP_waterVapour sinkP_waterVapour(pnom = 80e5, Tnom = 600, Res = 1) annotation(Placement(transformation(origin = {-46, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Blocks.Continuous.PI PI(k = 3, T = 100) annotation(Placement(transformation(extent = {{26, 24}, {44, 40}}, rotation = 0)));
      Modelica.Blocks.Math.Feedback feedback annotation(Placement(transformation(extent = {{-22, 22}, {-2, 42}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant step(k = 0.3) annotation(Placement(transformation(extent = {{-58, 22}, {-38, 42}}, rotation = 0)));
    equation
      connect(Evaporator.inlet, sourceW.flange) annotation(Line(points = {{16, -33.8}, {32, -33.8}, {32, -30}, {50, -30}}, color = {159, 159, 223}));
      connect(fixedHeatFlow.port, Evaporator.wall) annotation(Line(points = {{8, -8}, {8, -10.75}, {8, -10.75}, {8, -13.5}, {8, -36.9}, {8, -36.9}}, color = {191, 0, 0}));
      connect(sinkP_waterVapour.flange, Evaporator.outlet) annotation(Line(points = {{-36, -22}, {-20, -22}, {-20, -22}, {2, -22}}, color = {159, 159, 223}));
      connect(PI.y, sourceW.w) annotation(Line(points = {{44.9, 32}, {66.2, 32}, {66.2, -24}}, color = {0, 0, 127}));
      connect(feedback.y, PI.u) annotation(Line(points = {{-3, 32}, {24.2, 32}}, color = {0, 0, 127}));
      connect(Evaporator.voidFraction, feedback.u2) annotation(Line(points = {{0.4, -34}, {-12, -34}, {-12, 24}}, color = {0, 0, 127}));
      connect(step.y, feedback.u1) annotation(Line(points = {{-37, 32}, {-20, 32}}, color = {0, 0, 127}));
      annotation(Diagram(graphics));
    end WaterDrum_test;

    model WaterVapourFlow1d_test
      Components.Gas.SinkP sinkP(pnom = 80e5, Tnom = 100) annotation(Placement(transformation(origin = {-56, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Components.Water.Flow1DVapourPConstWall SuperHeater(V = 1, pstart = 90e5, Tstart = 1000, Tmstart = 1000) annotation(Placement(transformation(extent = {{6, -20}, {-14, 0}}, rotation = 0)));
      Components.Water.SourceW_waterVapour sourceW(wnom = 1, pnom = 90e5, Tnom = 800) annotation(Placement(transformation(extent = {{56, -20}, {36, 0}}, rotation = 0)));
      Components.Thermal.FixedHeatFlow fixedHeatFlow(Q_flow = 0) annotation(Placement(transformation(origin = {-4, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    equation
      connect(SuperHeater.inlet, sourceW.flange) annotation(Line(points = {{6, -10}, {36, -10}}, color = {159, 159, 223}));
      connect(sinkP.flange, SuperHeater.outlet) annotation(Line(points = {{-46, -10}, {-30, -10}, {-30, -10}, {-14, -10}}, color = {159, 159, 223}));
      connect(fixedHeatFlow.port, SuperHeater.wall) annotation(Line(points = {{-4, 10}, {-4, 6.25}, {-4, 6.25}, {-4, 2.5}, {-4, -5}, {-4, -5}}, color = {191, 0, 0}));
      annotation(Diagram(graphics));
    end WaterVapourFlow1d_test;

    model Turbine_test
      parameter Real fictitius_K = sourceW_waterVapour.wnom / sourceW_waterVapour.pnom;
      Components.Water.SourceW_waterVapour sourceW_waterVapour(pnom = 90e5, Tnom = 1000, wnom = 60) annotation(Placement(transformation(extent = {{-16, -38}, {4, -18}}, rotation = 0)));
      Components.Water.TurbineThermalPort turbineThermalPort(h_iso = 30e5, Gamma = 1, pstart = 90e5, k = fictitius_K) annotation(Placement(transformation(origin = {40, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Components.Thermal.FixedHeatFlow_DHT_N1 fixedHeatFlow_DHT_N1_1(Q_flow = 0) annotation(Placement(transformation(origin = {40, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    equation
      connect(fixedHeatFlow_DHT_N1_1.port, turbineThermalPort.dHT) annotation(Line(points = {{40, -14}, {40, -23.8}}, color = {255, 127, 0}));
      connect(sourceW_waterVapour.flange, turbineThermalPort.inlet) annotation(Line(points = {{4, -28}, {31, -28}}, color = {159, 159, 223}));
      annotation(Diagram(graphics));
    end Turbine_test;

    model Turbine_test_withThermalStress
      Components.Water.SourceW_waterVapour sourceW_waterVapour(pnom = 90e5, Tnom = 900, wnom = 70) annotation(Placement(transformation(extent = {{-60, -38}, {-40, -18}}, rotation = 0)));
      Components.Water.TurbineThermalPort turbineThermalPort(h_iso = 30e5, Gamma = 1, k = 1e-6) annotation(Placement(transformation(origin = {-8, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Components.Thermal.CylinderThermalStress_radial_noVector cylinderThermalStress(rint = 0.4, rext = 0.6, Tstartint = 800, Tstartext = 800) annotation(Placement(transformation(origin = {-8, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Components.Thermal.FixedHeatFlow_DHT_N1 fixedHeatFlow_DHT_N1_1(Q_flow = 0) annotation(Placement(transformation(origin = {-8, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
    equation
      connect(sourceW_waterVapour.flange, turbineThermalPort.inlet) annotation(Line(points = {{-40, -28}, {-17, -28}}, color = {159, 159, 223}));
      connect(turbineThermalPort.dHT, cylinderThermalStress.internalBoundary) annotation(Line(points = {{-8, -23.8}, {-8, -16.9}, {-8, -11}, {-8, -11}}, color = {255, 127, 0}));
      connect(fixedHeatFlow_DHT_N1_1.port, cylinderThermalStress.externalBoundary) annotation(Line(points = {{-8, 4}, {-8, 1.75}, {-8, 1.75}, {-8, -0.5}, {-8, -5}, {-8, -5}}, color = {255, 127, 0}));
      annotation(Diagram(graphics));
    end Turbine_test_withThermalStress;

    model braytonJoule_exhaustGases_test
      Components.Gas.BraytonJoule_exhaustGases braytonJoule_exhaustGases(T_min = 548, T_max = 843, w_min = 454, w_max = 585) annotation(Placement(transformation(extent = {{-26, -4}, {-6, 16}}, rotation = 0)));
      Components.Gas.SinkP sinkP annotation(Placement(transformation(extent = {{12, -4}, {32, 16}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 1, height = -1, offset = 1) annotation(Placement(transformation(extent = {{-70, -4}, {-50, 16}}, rotation = 0)));
    equation
      connect(braytonJoule_exhaustGases.flange, sinkP.flange) annotation(Line(points = {{-6, 6}, {-1.5, 6}, {-1.5, 6}, {3, 6}, {3, 6}, {12, 6}}, color = {159, 159, 223}));
      connect(ramp.y, braytonJoule_exhaustGases.load) annotation(Line(points = {{-49, 6}, {-43.25, 6}, {-43.25, 6}, {-37.5, 6}, {-37.5, 6}, {-26, 6}}, color = {0, 0, 127}));
      annotation(Diagram(graphics));
    end braytonJoule_exhaustGases_test;

    model SourceW_water_vapour_test
      Components.Water.SourceW_waterVapour sourceW_waterVapour(pnom = 90e5, Tnom = 1000, wnom = 60 * 10e4) annotation(Placement(transformation(extent = {{-60, -38}, {-40, -18}}, rotation = 0)));
      Components.Water.SinkP_waterVapour sinkP_waterVapour(pnom = 90e5, Tnom = 1000, Res = 0) annotation(Placement(transformation(extent = {{-14, -38}, {6, -18}}, rotation = 0)));
    equation
      connect(sourceW_waterVapour.flange, sinkP_waterVapour.flange) annotation(Line(points = {{-40, -28}, {-14, -28}}, color = {159, 159, 223}));
      annotation(Diagram(graphics));
    end SourceW_water_vapour_test;

    model Test_attemperation
      CombinedCycle.Optimization.Plants.CC_att_reference plant annotation(Placement(visible = true, transformation(origin = {-6, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      plant.load = 0.15;
      plant.TIT = 0;
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-6, Interval = 20));
    end Test_attemperation;

    model test_specific_heat
      CombinedCycle.Substances.WaterLiquid liquid;
      CombinedCycle.Substances.WaterVapour vapour;
      parameter Real p = 32.4e5;
      Real Tl;
      Real Tv;
      parameter Real Tminl = 32.9 + 273.15;
      parameter Real Tmaxl = 510;
      parameter Real Tminv = 512;
      parameter Real Tmaxv = 788;
      Real cpl;
      Real cpv;
    equation
      Tl = Tminl + (Tmaxl - Tminl) * time;
      Tv = Tminv + (Tmaxv - Tminv) * time;
      liquid.p = p;
      liquid.T = Tl;
      vapour.p = p;
      vapour.T = Tv;
      cpl = der(liquid.h) / (Tmaxl - Tminl);
      cpv = der(vapour.h) / (Tmaxv - Tminv);
    end test_specific_heat;
  end Circuits_Tests;

  package Optimization
    package Simulators
      model CC_Att_find_initial_state
        import Modelica.Utilities.Streams.*;
        String Path = "D:\\Pietro\\CCinit.txt";
        CombinedCycle.Optimization.Plants.CC_Att_WarmStartup plant(Attemperation = true) annotation(Placement(visible = true, transformation(origin = {-5, 9}, extent = {{-43, -43}, {43, 43}}, rotation = 0)));
      equation
        plant.load = 0.15;
        plant.TIT = 616.408;
        when terminal() then
          Modelica.Utilities.Files.remove(Path);
          print("  economizer.liquid[2].T = " + String(plant.economizer.liquid[2].T) + ";", Path);
          print("  economizer.liquid[3].T = " + String(plant.economizer.liquid[3].T) + ";", Path);
          print("  economizer.liquid[4].T = " + String(plant.economizer.liquid[4].T) + ";", Path);
          print("  economizer.liquid[5].T = " + String(plant.economizer.liquid[5].T) + ";", Path);
          print("  economizer_gas_side.gas[2].T = " + String(plant.economizer_gas_side.gas[2].T) + ";", Path);
          print("  economizer_gas_side.gas[3].T = " + String(plant.economizer_gas_side.gas[3].T) + ";", Path);
          print("  economizer_gas_side.gas[4].T = " + String(plant.economizer_gas_side.gas[4].T) + ";", Path);
          print("  economizer_gas_side.gas[5].T = " + String(plant.economizer_gas_side.gas[5].T) + ";", Path);
          print("  evaporator.alpha = " + String(plant.evaporator.alpha) + ";", Path);
          print("  evaporator.p = " + String(plant.evaporator.p) + ";", Path);
          print("  evaporator_gas_side.gas[2].T = " + String(plant.evaporator_gas_side.gas[2].T) + ";", Path);
          print("  evaporator_gas_side.gas[3].T = " + String(plant.evaporator_gas_side.gas[3].T) + ";", Path);
          print("  evaporator_gas_side.gas[4].T = " + String(plant.evaporator_gas_side.gas[4].T) + ";", Path);
          print("  PI_alpha.x = " + String(plant.PI_alpha.x) + "", Path);
          print("  superheater1.vapour[2].T = " + String(plant.superheater1.vapour[2].T) + ";", Path);
          print("  superheater1.vapour[3].T = " + String(plant.superheater1.vapour[3].T) + ";", Path);
          print("  superheater1.vapour[4].T = " + String(plant.superheater1.vapour[4].T) + ";", Path);
          print("  superheater1_gas_side.gas[2].T = " + String(plant.superheater1_gas_side.gas[2].T) + ";", Path);
          print("  superheater1_gas_side.gas[3].T = " + String(plant.superheater1_gas_side.gas[3].T) + ";", Path);
          print("  superheater1_gas_side.gas[4].T = " + String(plant.superheater1_gas_side.gas[4].T) + ";", Path);
          print("  superheater2.vapour[2].T = " + String(plant.superheater2.vapour[2].T) + ";", Path);
          print("  superheater2.vapour[3].T = " + String(plant.superheater2.vapour[3].T) + ";", Path);
          print("  superheater2.vapour[4].T = " + String(plant.superheater2.vapour[4].T) + ";", Path);
          print("  superheater2_gas_side.gas[2].T = " + String(plant.superheater2_gas_side.gas[2].T) + ";", Path);
          print("  superheater2_gas_side.gas[3].T = " + String(plant.superheater2_gas_side.gas[3].T) + ";", Path);
          print("  superheater2_gas_side.gas[4].T = " + String(plant.superheater2_gas_side.gas[4].T) + ";", Path);
          print("  turbineShaft.T[2:(turbineShaft.Nr-1)] = ones(turbineShaft.Nr-2) * " + String(plant.turbineShaft.T[2]) + ";", Path);
        end when;
        annotation(experiment(StopTime = 20000), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Documentation(info = "This model calculates the values of the state variables in the initial condition at steady state and prints the risult in a file specified by the String Path."));
      end CC_Att_find_initial_state;

      model StartupAttReference
        CombinedCycle.Optimization.Plants.CC_Att_WarmStartup plant annotation(Placement(visible = true, transformation(origin = {-5, 9}, extent = {{-43, -43}, {43, 43}}, rotation = 0)));
        parameter Real u0 = 0.15 "Initial value of control variable";
        parameter Real Tstart = 10000 "Start-up time";
        parameter Real N = 6;
        Real alpha;
        Real tau;
        Real k;
        //Real mu;
        Real lambda;
        Real cp;
        Real dh;
        Real dT;
        parameter Real sf = 1e2 "Smoothing factor (higher->steeper)";
        constant Real pi = Modelica.Constants.pi;
      equation
        dh = der(plant.superheater2.vapour[3].h);
        dT = der(plant.superheater2.vapour[3].T);
        cp = dh / dT;
        lambda = (atan(sf * (time - 4200)) + pi / 2) / pi;
        plant.load = u0 + (1 - u0) * (time / Tstart) / (1 + (time / Tstart) ^ N) ^ (1 / N);
        plant.TIT = 613 + 180 * (lambda + time / 4200 * (1 - lambda));
        alpha = (plant.superheater2_gas_side.Gi + plant.superheater2_gas_side.Gi + plant.superheater2_gas_side.Gi) / plant.superheater2.inlet.w / cp;
        tau = plant.superheater2.Cm / plant.superheater2.inlet.w / cp;
        k = exp(alpha) * Modelica.Constants.pi / 2 / tau;
        //mu = 2 * plant.superheater2.inlet.w * 2000 - (plant.superheater2_gas_side.node[1].G + plant.superheater2_gas_side.node[2].G + plant.superheater2_gas_side.node[3].G);
        annotation(experiment(StopTime = 10000), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Documentation(info = "This model replicates the reference simulation using in the optimization of the complete model of the Plant."));
      end StartupAttReference;
    end Simulators;



    package Plants
      model CC_Att_basic
        parameter Real w_nom = gasTurbine.w_max;
        parameter Integer Nv_eco(min = 1) = 4 "number of discrete volumes in the banks of the economizer";
        parameter Integer Nv_eva(min = 1) = 3 "number of discrete volumes in the banks of the evaporator";
        parameter Integer Nv_sh1(min = 1) = 3 "number of discrete volumes in the banks of the first superheater";
        parameter Integer Nv_sh2(min = 1) = 3 "number of discrete volumes in the banks of the second superheater";
        Modelica.Blocks.Interfaces.RealOutput sigma(nominal = 1e8) "External surface stress on rotor" annotation(Placement(transformation(extent = {{96, 50}, {116, 70}})));
        CombinedCycle.Components.Water.Flow1DVapourPConstNvolumes superheater1(N = Nv_sh1, Direction = false, Cm = 4e6 * 2 / 3, Tstart = 614, V = 5 * 2 / 3, pstart = 982455) annotation(Placement(visible = true, transformation(origin = {-12, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
        CombinedCycle.Components.Water.Flow1DVapourPConstNvolumes superheater2(N = Nv_sh2, Direction = false, Cm = 4e6 * 1 / 3, Tstart = 616, V = 5 * 1 / 3, pstart = 982455) annotation(Placement(visible = true, transformation(origin = {-60, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
        Modelica.Blocks.Interfaces.RealOutput p(nominal = 50e5) "External surface stress on rotor" annotation(Placement(visible = true, transformation(extent = {{96, -50}, {116, -30}}, rotation = 0), iconTransformation(extent = {{96, -70}, {116, -50}}, rotation = 0)));
        CombinedCycle.Components.Water.FlowJoin flowjoin annotation(Placement(visible = true, transformation(origin = {-36, 12}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        CombinedCycle.Components.Water.DrumNvolumes evaporator(pstart = 982455, alphastart = 0.5, V = 15, Cm = 200 * 70 / 4186 / 5, N = Nv_eva) annotation(Placement(visible = true, transformation(origin = {14, 2}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Components.Gas.BraytonJoule_exhaustGases gasTurbine(w_min = 115, wstart = 160, T_min = 529.15, T_max = 814, w_max = 200, load_min_w = 0.5, Tstart = 529.15) annotation(Placement(visible = true, transformation(extent = {{-80, -88}, {-60, -68}}, rotation = 0)));
        CombinedCycle.Components.Gas.FlowPConstNvolumes superheater1_gas_side(N = Nv_sh1, Gn = 200000 * 2 / 3, V = 4000 * 2 / 3, wn = w_nom, tau = 0.55, Tstart = 616) annotation(Placement(visible = true, transformation(origin = {-12, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CombinedCycle.Components.Gas.FlowPConstNvolumes superheater2_gas_side(N = Nv_sh2, Tstart = 617, V = 4000 * 1 / 3, Gn = 200000 * 1 / 3, wn = w_nom, tau = 0.55) annotation(Placement(visible = true, transformation(origin = {-42, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CombinedCycle.Components.Gas.FlowPConstNvolumes evaporator_gas_side(Gn = 700000, Tstart = 454, V = 4000, wn = w_nom, N = Nv_eva, tau = 0.50, pstart = 101325) annotation(Placement(visible = true, transformation(origin = {14, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Blocks.Continuous.PI PI_alpha(k = -1000, T = 0.01, y_start = 26.10) annotation(Placement(visible = true, transformation(extent = {{22, 74}, {40, 90}}, rotation = 0)));
        Blocks.Math.Feedback feedback_alpha annotation(Placement(visible = true, transformation(extent = {{-4, 72}, {16, 92}}, rotation = 0)));
        Blocks.Sources.Constant alpha_SP(k = 0.5) annotation(Placement(visible = true, transformation(extent = {{-30, 72}, {-10, 92}}, rotation = 0)));
        CombinedCycle.Blocks.Interfaces.RealInput TIT annotation(Placement(visible = true, transformation(origin = {-106, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
        CombinedCycle.Blocks.Interfaces.RealInput load annotation(Placement(visible = true, transformation(extent = {{-126, -98}, {-86, -58}}, rotation = 0), iconTransformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        Components.Water.TurbineThermalPort turbine(Gamma = 2000, k = 8.00e-006) annotation(Placement(visible = true, transformation(extent = {{-76, -10}, {-96, 10}}, rotation = 0)));
        Components.Thermal.CylinderThermalStressConstant turbineShaft(rint(displayUnit = "mm") = 0.04, rext(displayUnit = "mm") = 0.327, rho = 7750, c = 577, lambda = 43, linearExpansionCoefficient = 1.2e-5, youngModulus(displayUnit = "Pa") = 1.96e11, poissonRatio = 0.3, Nr = 8, Tstartint = 550, Tstartext = 550) annotation(Placement(visible = true, transformation(extent = {{-96, -8}, {-76, -28}}, rotation = 0)));
        CombinedCycle.Components.Water.Flow1DLiquid_inputw pseudo_valve(w(start = 0, fixed = true)) annotation(Placement(visible = true, transformation(origin = {-10, 28}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Components.Water.SourceW_waterLiquid_Input WaterLiquidSource(Tnom = 306.05, wstart = 26.10, pstart = 1002455) annotation(Placement(visible = true, transformation(extent = {{98, -14}, {78, 6}}, rotation = 0)));
        CombinedCycle.Components.Water.Flow1DLiquidPConstNvolumes economizer(pstart = 1002455, Cm = 5.89e6, V = 30, Direction = false, N = Nv_eco) annotation(Placement(visible = true, transformation(origin = {62, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
        Components.Gas.SinkP environment annotation(Placement(visible = true, transformation(extent = {{80, -88}, {100, -68}}, rotation = 0)));
        CombinedCycle.Components.Gas.FlowPConstNvolumes economizer_gas_side(Tstart = 422, V = 4000, Gn = 375000, wn = w_nom, N = Nv_eco, tau = 0.9) annotation(Placement(visible = true, transformation(origin = {62, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CombinedCycle.Components.Water.FlowSplit flowsplit annotation(Placement(visible = true, transformation(origin = {30, 24}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        CombinedCycle.Components.Water.Valve_DeltaPconst valve_economizer(deltaP = 2e5) annotation(Placement(visible = true, transformation(origin = {36, 4}, extent = {{-8, -8}, {8, 8}}, rotation = 90)));
      equation
        connect(valve_economizer.inlet, economizer.outlet) annotation(Line(points = {{36, -4}, {52, -4}, {52, -4}, {52, -4}}, color = {159, 159, 223}));
        connect(valve_economizer.outlet, flowsplit.inlet) annotation(Line(points = {{36, 12}, {36, 24}}, color = {159, 159, 223}));
        connect(evaporator.inlet, flowsplit.outlet2) annotation(Line(points = {{22, -4}, {25, -4}, {25, 20}, {24, 20}}, color = {159, 159, 223}));
        connect(flowsplit.outlet1, pseudo_valve.inlet) annotation(Line(points = {{24, 28}, {0, 28}}, color = {159, 159, 223}));
        connect(economizer_gas_side.wall, economizer.wall) annotation(Line(points = {{62, -73}, {62, -9}}, color = {255, 127, 0}));
        connect(evaporator_gas_side.outlet, economizer_gas_side.inlet) annotation(Line(points = {{24, -78}, {52, -78}}, color = {159, 159, 223}));
        connect(economizer_gas_side.outlet, environment.flange) annotation(Line(points = {{72, -78}, {80, -78}}, color = {159, 159, 223}));
        connect(WaterLiquidSource.flange, economizer.inlet) annotation(Line(points = {{78, -4}, {72, -4}}, color = {159, 159, 223}));
        connect(PI_alpha.y, WaterLiquidSource.w) annotation(Line(points = {{41, 82}, {94, 82}, {94, 2}}, color = {0, 0, 127}));
        connect(pseudo_valve.outlet, flowjoin.inlet1) annotation(Line(points = {{-20, 28}, {-30, 28}, {-30, 18}, {-30, 18}}, color = {159, 159, 223}));
        connect(turbine.shaft_wall, turbineShaft.externalBoundary) annotation(Line(points = {{-86, -4.2}, {-86, -15.2}}, color = {255, 127, 0}));
        connect(turbine.inlet, superheater2.outlet) annotation(Line(points = {{-79, 8}, {-70, 8}}, color = {159, 159, 223}));
        connect(alpha_SP.y, feedback_alpha.u1) annotation(Line(points = {{-9, 82}, {-2, 82}}, color = {0, 0, 127}));
        connect(feedback_alpha.y, PI_alpha.u) annotation(Line(points = {{15, 82}, {17.6, 82}, {20.2, 82}}, color = {0, 0, 127}));
        connect(evaporator.voidFraction, feedback_alpha.u2) annotation(Line(points = {{6, -4}, {6, 74}}, color = {0, 0, 127}));
        connect(superheater1_gas_side.outlet, evaporator_gas_side.inlet) annotation(Line(points = {{-2, -78}, {4, -78}}, color = {159, 159, 223}));
        connect(evaporator_gas_side.wall, evaporator.wall) annotation(Line(points = {{14, -73}, {14, -7}}, color = {255, 127, 0}));
        connect(gasTurbine.flange, superheater2_gas_side.inlet) annotation(Line(points = {{-60, -78}, {-52, -78}}, color = {159, 159, 223}));
        connect(superheater2_gas_side.wall, superheater2.wall) annotation(Line(points = {{-42, -73}, {-42, -34}, {-60, -34}, {-60, 3}}, color = {255, 238, 44}));
        connect(superheater2_gas_side.outlet, superheater1_gas_side.inlet) annotation(Line(points = {{-32, -78}, {-22, -78}}, color = {159, 159, 223}));
        connect(superheater1.wall, superheater1_gas_side.wall) annotation(Line(points = {{-12, 3}, {-12, -73}}, color = {255, 238, 44}));
        connect(gasTurbine.load, load) annotation(Line(points = {{-80, -78}, {-106, -78}}, color = {0, 0, 127}));
        connect(superheater1.inlet, evaporator.outlet) annotation(Line(points = {{-2, 8}, {8, 8}}, color = {159, 159, 223}));
        connect(flowjoin.inlet2, superheater1.outlet) annotation(Line(points = {{-30, 8}, {-22, 8}}, color = {159, 159, 223}));
        connect(flowjoin.outlet, superheater2.inlet) annotation(Line(points = {{-42, 12}, {-50, 12}, {-50, 8}}, color = {159, 159, 223}));
        // the minus sign is used to obtain a positive value
        // the factor 2 is the geometric factor for thermal
        // stresses as in Casella, Pretolani 2006.
        sigma = -2 * turbineShaft.extThermalStress;
        p = evaporator.p;
        annotation(Diagram(graphics), experiment(StopTime = 1000), experimentSetupOutput, Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2}), graphics = {Rectangle(lineColor = {0, 0, 255}, extent = {{-100, 100}, {100, -100}}), Text(lineColor = {0, 0, 255}, extent = {{-136, 78}, {138, -74}}, textString = "P"), Text(origin = {-105, 34}, extent = {{-15, 6}, {25, -14}}, textString = "Attemperation"), Text(origin = {-102, -28}, extent = {{-18, 8}, {22, -12}}, textString = "Load")}), Documentation(info = "<HTML>
      <p>Partial model of the plant with the components common to both configurations.
      </HTML>"));
      end CC_Att_basic;
    
      model CC_Att
        extends CombinedCycle.Optimization.Plants.CC_Att_basic;
        CombinedCycle.Blocks.Pseudo_controller pseudo_controller(Attemperation = Attemperation) annotation(Placement(visible = true, transformation(origin = {-18, 48}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
        CombinedCycle.Blocks.Math.Feedback feedback_TIT annotation(Placement(visible = true, transformation(origin = {-72, 50}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
        parameter Boolean Attemperation = true "Select false to deactivate the attemperation";
        CombinedCycle.Blocks.Continuous.I I(k = 0.08, y_start = 616.254) annotation(Placement(visible = true, transformation(origin = {-46, 50}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
        Real w_att;
      equation
        w_att = pseudo_valve.w;
        connect(I.y, pseudo_controller.T0) annotation(Line(points = {{-39, 50}, {-32.5, 50}, {-32.5, 48}, {-26, 48}}, color = {0, 0, 127}));
        connect(feedback_TIT.y, I.u) annotation(Line(points = {{-61, 50}, {-53, 50}}, color = {0, 0, 127}));
        connect(TIT, feedback_TIT.u1) annotation(Line(points = {{-106, 50}, {-82, 50}}, color = {0, 0, 127}));
        connect(turbine.TIT, feedback_TIT.u2) annotation(Line(points = {{-82, 9}, {-82, 24.5}, {-72, 24.5}, {-72, 40}}, color = {0, 0, 127}));
        connect(flowjoin.Tout, pseudo_controller.Tm) annotation(Line(points = {{-39, 15}, {-38, 15}, {-38, 40}, {-18, 40}}, color = {0, 0, 127}));
        connect(pseudo_controller.w, pseudo_valve.w) annotation(Line(points = {{-10, 48}, {-10, 34}}, color = {0, 0, 127}));
      annotation(
          Documentation(info = "Full model of the plant.
It has two control variables, the load of the gas turbine and the set point of the temparature of the steam at the inlet of the turbine controlled by the attemperation."));end CC_Att;

      model CC_att_reference
        extends CombinedCycle.Optimization.Plants.CC_Att_basic;
      equation
        connect(TIT, pseudo_valve.w) annotation(Line(points = {{-106, 50}, {-10, 50}, {-10, 34}, {-10, 34}}, color = {0, 0, 127}));
      annotation(
          Documentation(info = "Reference model for the plant.
In this configuration, the attemperation is not used for the control of the plant."));end CC_att_reference;

      model CC_Att_WarmStartup
        extends CombinedCycle.Optimization.Plants.CC_Att;
      initial equation
        economizer.liquid[2].T = 450.652;
        economizer.liquid[3].T = 452.023;
        economizer.liquid[4].T = 452.417;
        economizer.liquid[5].T = 452.53;
        economizer_gas_side.gas[2].T = 452.527;
        economizer_gas_side.gas[3].T = 452.405;
        economizer_gas_side.gas[4].T = 451.981;
        economizer_gas_side.gas[5].T = 422.734;
        evaporator.alpha = 0.5;
        evaporator.p = 982863;
        evaporator_gas_side.gas[2].T = 470.21;
        evaporator_gas_side.gas[3].T = 454.582;
        evaporator_gas_side.gas[4].T = 452.562;
        PI_alpha.x = -0.00786291;
        superheater1.vapour[2].T = 595.822;
        superheater1.vapour[3].T = 614.56;
        superheater1.vapour[4].T = 616.254;
        superheater1_gas_side.gas[2].T = 616.148;
        superheater1_gas_side.gas[3].T = 613.409;
        superheater1_gas_side.gas[4].T = 591.127;
        superheater2.vapour[2].T = 616.351;
        superheater2.vapour[3].T = 616.391;
        superheater2.vapour[4].T = 616.408;
        superheater2_gas_side.gas[2].T = 616.415;
        superheater2_gas_side.gas[3].T = 616.41;
        superheater2_gas_side.gas[4].T = 616.395;
        turbineShaft.T[2:turbineShaft.Nr - 1] = ones(turbineShaft.Nr - 2) * 616.408;
        //I.y = 616.254;
        annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
          Documentation(info = "Full model of the plant with initial equation for the warm start up."));
      end CC_Att_WarmStartup;

      model CC_Att_reference_WarmStartup
        extends CombinedCycle.Optimization.Plants.CC_att_reference;
      initial equation
        economizer.liquid[2].T = 450.652;
        economizer.liquid[3].T = 452.023;
        economizer.liquid[4].T = 452.417;
        economizer.liquid[5].T = 452.53;
        economizer_gas_side.gas[2].T = 452.527;
        economizer_gas_side.gas[3].T = 452.405;
        economizer_gas_side.gas[4].T = 451.981;
        economizer_gas_side.gas[5].T = 422.734;
        evaporator.alpha = 0.5;
        evaporator.p = 982863;
        evaporator_gas_side.gas[2].T = 470.21;
        evaporator_gas_side.gas[3].T = 454.582;
        evaporator_gas_side.gas[4].T = 452.562;
        PI_alpha.x = -0.00786291;
        superheater1.vapour[2].T = 595.822;
        superheater1.vapour[3].T = 614.56;
        superheater1.vapour[4].T = 616.254;
        superheater1_gas_side.gas[2].T = 616.148;
        superheater1_gas_side.gas[3].T = 613.409;
        superheater1_gas_side.gas[4].T = 591.127;
        superheater2.vapour[2].T = 616.351;
        superheater2.vapour[3].T = 616.391;
        superheater2.vapour[4].T = 616.408;
        superheater2_gas_side.gas[2].T = 616.415;
        superheater2_gas_side.gas[3].T = 616.41;
        superheater2_gas_side.gas[4].T = 616.395;
        turbineShaft.T[2:turbineShaft.Nr - 1] = ones(turbineShaft.Nr - 2) * 616.408;
        annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})),
          Documentation(info = "Reference model of the plant with initial equation for the warm start up."));
      end CC_Att_reference_WarmStartup;
    end Plants;


  end Optimization;

  package Types
    type Temperature = Modelica.SIunits.Temperature(nominal = 500);
    type Pressure = Modelica.SIunits.Pressure(nominal = 1e6);
    type GasDensity = Modelica.SIunits.Density(nominal = 30);
    type VapourDensity = Modelica.SIunits.Density(nominal = 30);
    type LiquidDensity = Modelica.SIunits.Density(nominal = 1000);
    type MassFlowRate = Modelica.SIunits.MassFlowRate(nominal = 10);
    type HeatFlowRate = Modelica.SIunits.HeatFlowRate(nominal = 1e6);
    type SpecificEnthalpy = Modelica.SIunits.SpecificEnthalpy;
    type SpecificEnthalpyGas = Modelica.SIunits.SpecificEnthalpy(nominal = 1e6);
    type SpecificEnergyGas = Modelica.SIunits.SpecificEnergy(nominal = 1e6);
    type SpecificEnthalpyLiquid = Modelica.SIunits.SpecificEnthalpy(nominal = 1e6);
    type SpecificEnergyLiquid = Modelica.SIunits.SpecificEnergy(nominal = 1e6);
    type SpecificHeatCapacity = Modelica.SIunits.SpecificHeatCapacity(nominal = 1000);
    type HeatFlux = Real(final quantity = "HeatFlux", final unit = "W/m2", nominal = 1e4);
    type Stress = Modelica.SIunits.NormalStress(nominal = 1e7);
    type VolumetricEnergyLiquid = Real(final quantity = "VolumetricEnergyLiquid", final unit = "J/(m3.Pa)");
    type VolumetricEnergyGas = Real(final quantity = "VolumetricEnergyGas", final unit = "J/(m3.Pa)");
    type DerSaturationTemperatureByPressure = Real(final quantity = "DerSaturationTemperatureByPressure", final unit = "K/Pa");
    type DerDensityByPressureLiquid = Real(final quantity = "DerDensityByPressureLiquid", final unit = "kg/(m3.Pa)");
    type DerDensityByPressureGas = Real(final quantity = "DerDensityByPressureGas", final unit = "kg/(m3.Pa)");
    type DerVolumetricEnergyByPressureLiquid = Real(final quantity = "DerVolumetricEnergyByPressureLiquid", final unit = "J/(m3.Pa)");
    type DerVolumetricEnergyByPressureGas = Real(final quantity = "DerVolumetricEnergyByPressureGas", final unit = "J/(m3.Pa)");
    type ThermalConductance = Modelica.SIunits.ThermalConductance(nominal = 1e6);
    type VapourMass = Modelica.SIunits.Mass(nominal = 30 * 7.5);
    type LiquidMass = Modelica.SIunits.Mass(nominal = 700 * 7.5);
  end Types;

  package test
    model economizer
      parameter Real w_nom = 200.0;
      parameter Integer Nv_eco(min = 1) = 4 "number of discrete volumes in the banks of the economizer";
      parameter Integer Nv_eva(min = 1) = 3 "number of discrete volumes in the banks of the evaporator";
      parameter Integer Nv_sh1(min = 1) = 3 "number of discrete volumes in the banks of the first superheater";
      parameter Integer Nv_sh2(min = 1) = 3 "number of discrete volumes in the banks of the second superheater";
      Components.Water.SourceW_waterLiquid_Input WaterLiquidSource(Tnom = 306.05, wstart = 26.10, pstart = 1002455) annotation(Placement(visible = true, transformation(extent = {{98, -14}, {78, 6}}, rotation = 0)));
      CombinedCycle.Components.Water.Flow1DLiquidPConstNvolumes economizer(pstart = 1002455, Cm = 5.89e6, V = 30, Direction = false, N = Nv_eco) annotation(Placement(visible = true, transformation(origin = {62, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      CombinedCycle.Components.Gas.FlowPConstNvolumes economizer_gas_side(Tstart = 422, V = 4000, Gn = 375000, wn = w_nom, N = Nv_eco, tau = 0.9) annotation(Placement(visible = true, transformation(origin = {62, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Components.Gas.SinkP environment annotation(Placement(visible = true, transformation(extent = {{80, -88}, {100, -68}}, rotation = 0)));
      Components.Gas.SinkP sinkP1(Res = 0, pnom = 1.18286e+06) annotation(Placement(visible = true, transformation(origin = {22, -4}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Components.Gas.SourceW sourceW1(G = 0, wnom = 114.459, Tnom = 452.5618) annotation(Placement(visible = true, transformation(origin = {30, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      WaterLiquidSource.w = 7.86291;
      connect(sourceW1.flange, economizer_gas_side.inlet) annotation(Line(points = {{40, -78}, {52, -78}, {52, -78}, {52, -78}}, color = {159, 159, 223}));
      connect(sinkP1.flange, economizer.outlet) annotation(Line(points = {{32, -4}, {52, -4}, {52, -4}, {52, -4}}, color = {159, 159, 223}));
      connect(economizer_gas_side.outlet, environment.flange) annotation(Line(points = {{72, -78}, {80, -78}, {80, -78}, {80, -78}}, color = {159, 159, 223}));
      connect(economizer_gas_side.wall, economizer.wall) annotation(Line(points = {{62, -73}, {62, -73}, {62, -8}, {62, -8}}, color = {255, 238, 44}));
      connect(WaterLiquidSource.flange, economizer.inlet) annotation(Line(points = {{78, -4}, {72, -4}, {72, -4}, {72, -4}}, color = {159, 159, 223}));
      /*initial equation
                                                                                                                    economizer.liquid[2].T = 350.652;
                                                                                                                    economizer.liquid[3].T = 352.023;
                                                                                                                    economizer.liquid[4].T = 352.417;
                                                                                                                    economizer.liquid[5].T = 352.53;
                                                                                                                    economizer_gas_side.gas[2].T = 452.527;
                                                                                                                    economizer_gas_side.gas[3].T = 452.405;
                                                                                                                    economizer_gas_side.gas[4].T = 451.981;
                                                                                                                    economizer_gas_side.gas[5].T = 422.734;*/
      annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
    end economizer;

    model test
      parameter Real w_nom = 200.0;
      parameter Integer Nv_eco(min = 1) = 4 "number of discrete volumes in the banks of the economizer";
      parameter Integer Nv_eva(min = 1) = 3 "number of discrete volumes in the banks of the evaporator";
      parameter Integer Nv_sh1(min = 1) = 3 "number of discrete volumes in the banks of the first superheater";
      parameter Integer Nv_sh2(min = 1) = 3 "number of discrete volumes in the banks of the second superheater";
      Components.Water.SourceW_waterLiquid_Input WaterLiquidSource(Tnom = 306.05, wstart = 26.10, pstart = 1002455) annotation(Placement(visible = true, transformation(extent = {{98, -14}, {78, 6}}, rotation = 0)));
      CombinedCycle.Components.Water.Flow1DLiquidPConstNvolumes economizer(pstart = 1.18286e+06, Cm = 5.89e6, V = 30, Direction = false, N = Nv_eco) annotation(Placement(visible = true, transformation(origin = {62, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Components.Gas.SinkP sinkP1(Res = 0, pnom = 1.18286e+06) annotation(Placement(visible = true, transformation(origin = {22, -4}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Components.Thermal.Wall_fixedT wall_fixedT1(T = 400, N = Nv_eco, G = 1000) annotation(Placement(visible = true, transformation(origin = {60, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(wall_fixedT1.wall, economizer.wall) annotation(Line(points = {{60, -44}, {60, -44}, {60, -10}, {60, -10}}, color = {255, 238, 44}));
      WaterLiquidSource.w = 7.86291;
      connect(sinkP1.flange, economizer.outlet) annotation(Line(points = {{32, -4}, {52, -4}, {52, -4}, {52, -4}}, color = {159, 159, 223}));
      connect(WaterLiquidSource.flange, economizer.inlet) annotation(Line(points = {{78, -4}, {72, -4}, {72, -4}, {72, -4}}, color = {159, 159, 223}));
    initial equation
      economizer.liquid[2].T = 306.05;
      economizer.liquid[3].T = 306.05;
      economizer.liquid[4].T = 306.05;
      economizer.liquid[5].T = 306.05;
    end test;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end test;
  annotation(Modelica(version = "2.2.1"), uses(Modelica(version = "3.1")), version = "2", conversion(from(version = "", script = "ConvertFromCombinedCycle_.mos"), noneFromVersion = "1"));
end CombinedCycle;
