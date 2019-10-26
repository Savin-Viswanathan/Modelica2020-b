package OceanEngineering
  extends Modelica.Icons.Package;

  package Components
    package Waves
      package RegularWave
        model Regular_Airy_Wave
          extends Modelica.Blocks.Icons.Block;
          OceanEngineering.Connectors.WaveDataConnector wdoc(n_omega_i = 1) annotation(
            Placement(transformation(extent = {{80, -20}, {120, 20}})));
          constant Real pi = Modelica.Constants.pi "Value of pi";
          constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
          parameter Modelica.SIunits.Length d = 100 "Water depth";
          parameter Modelica.SIunits.Density rho_w = 1025 "Density of Water";
          parameter Modelica.SIunits.Length Hr = 1 "Wave Height";
          parameter Real Tr = 7 "Wave Period";
          parameter Real Trmp = 20 "Wave ramp time";
          parameter Real Tdel = 0 "Wave delay";
          parameter Integer n_omega_i = 1 "Number of frequency components";
          Real T[size(omega, 1)] "Wave period";
          Modelica.SIunits.AngularFrequency omega[n_omega_i] "Frequency components of time series";
          Modelica.SIunits.Length zeta0i[n_omega_i] "Amplitude of wave component";
          Modelica.SIunits.Length SSE_X0 "Sea surface elevation at x=0 - to plot generated waves";
          Real epsilon[n_omega_i] "Random phase of wave components";
          Real k[n_omega_i](each unit = "m^-1") = OceanEngineering.Functions.waveNumberIterator(d, omega) "Wave number of component waves";
        equation
          for i in 1:n_omega_i loop
            omega[i] = 2 * pi / Tr;
            epsilon[i] = 0;
            if time < Tdel then
              zeta0i[i] = 0;
            elseif time < Tdel + Trmp then
              zeta0i[i] = sin(pi / 2 * (time - Tdel) / Trmp) * Hr / 2;
            else
              zeta0i[i] = Hr / 2;
            end if;
            T[i] = Tr;
          end for;
          SSE_X0 = sum(zeta0i .* cos(omega * time - 2 * pi * epsilon));
          wdoc.d = d;
          wdoc.rho_w = rho_w;
          wdoc.omega = omega;
          wdoc.T = T;
          wdoc.k = k;
          wdoc.epsilon = epsilon;
          wdoc.zeta0i = zeta0i;
          wdoc.SSE_X0 = SSE_X0;
          annotation(
            Icon(graphics = {Text(origin = {-54, 65}, extent = {{-40, 17}, {40, -17}}, textString = "H,T"), Rectangle(origin = {100, 0}, fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}}), Line(origin = {-15, 0}, points = {{-81, 0}, {81, 0}}), Line(origin = {-14.97, 1.01}, points = {{-81.0281, -1.00867}, {-33.0281, 42.9913}, {36.9719, -43.0087}, {80.9719, -1.00867}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
            experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.5));
        end Regular_Airy_Wave;


      end RegularWave;

      package IrregularWave
        model IRW_PM_RDFCWI
          extends Modelica.Blocks.Icons.Block;
          OceanEngineering.Connectors.WaveDataConnector wdoc(n_omega_i = 100) annotation(
            Placement(transformation(extent = {{80, -20}, {120, 20}})));
          constant Real pi = Modelica.Constants.pi "Value of pi";
          constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
          parameter Modelica.SIunits.Length d = 100 "Water depth";
          parameter Modelica.SIunits.Density rho_w = 1025 "Density of Water";
          parameter Modelica.SIunits.Length Hs = 1 "Significant Wave Height";
          parameter Modelica.SIunits.AngularFrequency omega_min = 0.03141 "Lowest frequency component/frequency interval";
          parameter Modelica.SIunits.AngularFrequency omega_max = 3.141 "Highest frequency component";
          parameter Integer n_omega_i = 100 "Number of frequency components";
          parameter Integer localSeed = 614657 "Local seed for random number generator";
          parameter Integer globalSeed = 30020 "Global seed for random number generator";
          parameter Real rnd_shft[n_omega_i] = OceanEngineering.Functions.randomNumberGenerator(localSeed, globalSeed, n_omega_i);
          parameter Integer localSeed1 = 614757 "Local seed for random number generator";
          parameter Integer globalSeed1 = 40020 "Global seed for random number generator";
          parameter Real epsilon[n_omega_i] = OceanEngineering.Functions.randomNumberGenerator(localSeed1, globalSeed1, n_omega_i);
          parameter Integer SSE_X0on = 1 "Flag for generating SSE profile";
          parameter Real Trmp = 200 "Interval for ramping up of waves during start phase";
          parameter Real omega[n_omega_i] = OceanEngineering.Functions.frequencySelector(omega_min, omega_max, rnd_shft);
          parameter Real S[n_omega_i] = OceanEngineering.Functions.spectrumGenerator_PM(Hs, omega);
          parameter Modelica.SIunits.Length zeta0i[n_omega_i] = sqrt(2 * S * omega_min) "Amplitude of wave component";
          parameter Real T[n_omega_i] = 2 * pi ./ omega "Wave period";
          parameter Real k[n_omega_i] = OceanEngineering.Functions.waveNumberIterator(d, omega) "Wave number";
          Modelica.SIunits.Length SSE_X0 "Sea surface elevation at x=0 - to plot generated waves";
          Real zeta0i_rmp[n_omega_i] "Ramp value of zeta0i";
        equation
          for i in 1:n_omega_i loop
            if time < Trmp then
              zeta0i_rmp[i] = sin(pi / 2 * time / Trmp) * zeta0i[i];
            else
              zeta0i_rmp[i] = zeta0i[i];
            end if;
          end for;
          if SSE_X0on == 1 then
            SSE_X0 = sum(zeta0i_rmp .* cos(omega * time - 2 * pi * epsilon));
          else
            SSE_X0 = 0;
          end if;
          wdoc.d = d;
          wdoc.rho_w = rho_w;
          wdoc.omega = omega;
          wdoc.T = T;
          wdoc.k = k;
          wdoc.epsilon = epsilon;
          wdoc.zeta0i = zeta0i_rmp;
          wdoc.SSE_X0 = SSE_X0;
          annotation(
            Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
            experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.5));
        end IRW_PM_RDFCWI;
      end IrregularWave;
    end Waves;

    package CurrentProfiles
      block CurrentProfile_4pt
        extends Modelica.Blocks.Icons.Block;
        OceanEngineering.Connectors.CurrentDataConnector cdc annotation(
          Placement(transformation(extent = {{80, -20}, {120, 20}})));
        constant Real pi = Modelica.Constants.pi "Value of pi";
        parameter Real zcg[:] = {-100, -80, -10, 0};
        parameter Real Uf[:] = {0, 0.5, 1, 1.5};
        parameter Real Trmp = 10;
        Real Ucg[size(Uf, 1)];
      equation
        for i in 1:size(Uf, 1) loop
          if time < Trmp then
            Ucg[i] = sin(pi / 2 * time / Trmp) * Uf[i];
          else
            Ucg[i] = Uf[i];
          end if;
        end for;
        cdc.zcg = zcg;
        cdc.Ucg = Ucg;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -30}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 70}, {100, -70}}), Line(origin = {-80, -29}, points = {{0, 69}, {0, -69}}), Line(origin = {-50.0068, -29.0068}, points = {{62.0068, 69.0068}, {10.0068, 17.0068}, {-9.99322, -58.9932}, {-29.9932, -68.9932}}), Line(origin = {-53, 16.0755}, points = {{-27, 0}, {41, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-60, -14}, points = {{-20, 0}, {20, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-65, -48}, points = {{-15, 0}, {15, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-70, -86}, points = {{-10, 0}, {10, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-47.6886, 40.3302}, points = {{-32, 0}, {60, 0}}, arrow = {Arrow.None, Arrow.Filled}), Rectangle(origin = {100, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}})}, coordinateSystem(initialScale = 0.1)));
      end CurrentProfile_4pt;
    end CurrentProfiles;

    package Floaters
      model CylindricalBuoy
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_a TopHook(s = z2) "Hook for spring" annotation(
          Placement(transformation(extent = {{-20, 120}, {20, 80}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a Fairlead[2](s = {x, z3}) annotation(
          Placement(transformation(extent = {{-20, -120}, {20, -80}})));
        OceanEngineering.Connectors.EnvironmentBuoyDataConnector ebdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{-80, -20}, {-120, 20}}, rotation = 0), iconTransformation(extent = {{-118, -20}, {-78, 20}}, rotation = 0)));
        constant Real g = Modelica.Constants.g_n "Value of acceleration due to gravity";
        constant Real pi = Modelica.Constants.pi "Value of pi";
        parameter Integer nOmega = 1 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        parameter Real d = 100;
        parameter Real rho_w = 1025;
        parameter Modelica.SIunits.Length r = 0.6 "Radius of the buoy";
        parameter Modelica.SIunits.Length h = 2 "Height of the buoy";
        parameter Modelica.SIunits.Mass m_s = 350 "Structural mass of the buoy";
        parameter Modelica.SIunits.Mass m_b = 500 "Ballast mass";
        parameter Modelica.SIunits.Length pz = 2 "Depth of panel to calculate current force & Morison Force";
        parameter Modelica.SIunits.Length zKG = 1.25 "KG of buoy";
        parameter Real Cma_x = 1 "Added mass coefficient";
        parameter Real Cma_z = 1 "Added mass coefficient";
        parameter Real Cd_x = 1 "Drag coefficient of buoy";
        parameter Real Cd_z = 1 "Drag coefficient of buoy";
        parameter Real K_x(unit = "N m^-1") = 0 "Coefficient of Stiffness";
        parameter Real K_z(unit = "N m^-1") = rho_w * g * pi * r ^ 2 "Coefficient of Stiffness";
        parameter Real C_x(unit = "kg/s") = 0 "Damping in x";
        parameter Real C_z(unit = "kg/s") = 3100 "Damping in z";
        parameter Modelica.SIunits.Mass M = m_s + m_b "Dry mass of the buoy";
        parameter Modelica.SIunits.Area Awp = pi * r ^ 2 "Waterplane area of buoy";
        parameter Integer mor_flg = 1 " '1' if moored, '0' if free floating";
        parameter Modelica.SIunits.Length z_s = M / (Awp * rho_w) "Static draught";
        parameter Modelica.SIunits.Length z_fb = zKG - z_s "VCG of body wrt SWL";
        parameter Modelica.SIunits.Mass spm_chain = 14 "Specific mass of the mooring line";
        parameter Modelica.SIunits.Density rho_mat_chain = 7800 "Density of mooring line material";
        parameter Modelica.SIunits.Mass spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Modelica.SIunits.Length z_m = z_s+mor_flg*(spm_chain_sub * (d-z_s) * g / K_z) "Approximate change in static draught when mooring chain is attached";
        parameter Modelica.SIunits.Length z_cb = -z_m / 2 "Static centre of buoyancy";
        Modelica.SIunits.Length SSE_X "Sea surface elevation corresponding to buoy CG X coordinate";
        Modelica.SIunits.Acceleration a_w "Water particle acceleration in z direction";
        Modelica.SIunits.Length x "x co-ordinate of VCG";
        Modelica.SIunits.Length z "Displacement of body VCG from z_fb";
        Modelica.SIunits.Velocity v_x "Velocity of body VCG in horizontal direction";
        Modelica.SIunits.Acceleration a_x "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Velocity v_z "Velocity of body VCG in vertical direction";
        Modelica.SIunits.Acceleration a_z "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Length z1 "Position of body VCG";
        Modelica.SIunits.Length z2 "Position of deck";
        Modelica.SIunits.Length z3 "Position of keel";
        Modelica.SIunits.Force MF_x "Morison force due to current and waves";
      initial equation
        z3=-z_m;
        der(z) = 0;
        x = ebdc.xinit;
        der(x) = 0;
      equation
        SSE_X = sum(ebdc.zeta0i .* cos(ebdc.k * x - ebdc.omega * time - 2 * pi * ebdc.epsilon));
        a_w = OceanEngineering.Functions.wave_awCalculator(time, d, ebdc.omega, ebdc.T, ebdc.k, ebdc.zeta0i, ebdc.epsilon, x, z_cb, SSE_X);
        v_z = der(z);
        a_z = der(v_z);
        M * a_z + C_z * (SSE_X - z3) * v_z + K_z * z = rho_w * g * Awp * SSE_X + Cma_x * M * a_w + TopHook.f + Fairlead[2].f;
        z1 = z_fb + z;
        z2 = z3 + h;
        z3 = -z_s + z;
        MF_x = OceanEngineering.Functions.morisonForceCydlBuoy(time, d, ebdc.omega, ebdc.T, ebdc.k, x, v_x, a_x, ebdc.epsilon, ebdc.zeta0i, rho_w, SSE_X, r, h, Cd_x, Cma_x, z3, pz, ebdc.zcg, ebdc.Ucg);
        v_x = der(x);
        a_x = der(v_x);
        M * a_x + C_x * M * (SSE_X - z3) * v_x + K_x * x = Fairlead[1].f + MF_x;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -50}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 50}, {100, -50}}), Rectangle(origin = {0, 3}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 53}, {28, -53}}), Ellipse(origin = {0, 57}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360), Ellipse(origin = {0, -51}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-05, Interval = 0.5));
      end CylindricalBuoy;










    end Floaters;

    package Moorings
      model CatenaryMooring_Mf0
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_b shackle[2] annotation(
          Placement(transformation(extent = {{-20, 80}, {20, 120}})));
        OceanEngineering.Connectors.EnvironmentMooringDataConnector emdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{-80, -20}, {-120, 20}}, rotation = 0), iconTransformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        parameter Integer nOmega = 100 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        constant Real g = Modelica.Constants.g_n;
        constant Real pi = Modelica.Constants.pi;
        parameter Real rho_w = 1025;
        parameter Real d = 100;
        parameter Real lnk_l = 3;
        parameter Integer n_lnk = 50;
        parameter Real L_chain = lnk_l * n_lnk;
        parameter Real d_chain = 0.022;
        parameter Real rho_mat_chain = 7800;
        parameter Real spm_chain = 10;
        parameter Real spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Real Xmax = sqrt(L_chain ^ 2 - d ^ 2);
        parameter Real X0 = L_chain - d;
        parameter Real xmax = Xmax - 10;
        parameter Real x[:] = 0:0.1:xmax;
        parameter Real Th[:] = OceanEngineering.Functions.catThIterator(L_chain, x, d, spm_chain_sub);
        parameter Real X[:] = OceanEngineering.Functions.catXIterator(L_chain, d, spm_chain_sub, Th, x);
        parameter Real Cdn = 0.75;
        parameter Real Cdt = 0;
        parameter Real Cman = 1;
        parameter Real Cmat = 0.0;
        parameter Real Cmn = 1 + Cman;
        parameter Real Cmt = 1 + Cmat;
        parameter Real D = 1.8 * d_chain;
        parameter Real Ai = rho_w * pi * D ^ 2 / 4;
        parameter Real Ad = 0.5 * rho_w * D;
        parameter Real Trmp = 50;
        Real x_lnk_cat[n_lnk + 1];
        Real x_lnk_plot_o;
        Real z_lnk_cat[n_lnk + 1];
        Real x_lnk_plot[n_lnk + 1];
        Real z_lnk_plot[n_lnk + 1];
        Real fd0;
        Real a;
        Real z_cat;
        Real x_cat;
        Real psi_cat;
        Real s_cat;
      equation
        if shackle[1].s <= X0 - shackle[2].s then
          fd0 = 0;
          a = 0;
          z_cat = d + shackle[2].s;
          psi_cat = pi / 2;
          x_cat = X0 - shackle[2].s;
          s_cat = d + shackle[2].s;
          x_lnk_cat[1] = 0;
          x_lnk_plot_o = shackle[1].s;
          x_lnk_plot[1] = shackle[1].s;
          z_lnk_cat[1] = 0;
          z_lnk_plot[1] = shackle[2].s;
          for i in 2:n_lnk loop
            if z_lnk_plot[i - 1] > (-(d - lnk_l)) then
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = shackle[1].s;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = z_lnk_plot[i - 1] - lnk_l;
            elseif z_lnk_plot[i - 1] > (-d) then
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - sqrt(lnk_l ^ 2 - (d + z_lnk_plot[i - 1]) ^ 2);
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            else
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            end if;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          shackle[1].f = -10;
          shackle[2].f = s_cat * spm_chain_sub * g;
        elseif shackle[1].s < sqrt(L_chain ^ 2 - d ^ 2) then
          fd0 = OceanEngineering.Functions.linearInterpolatorSV(X, Th, shackle[1].s);
          a = fd0 / (spm_chain_sub * g);
          z_cat = a + d + shackle[2].s;
          psi_cat = acos(a / z_cat);
          x_cat = fd0 / (spm_chain_sub * g) * Modelica.Math.acosh(1 + spm_chain_sub * g * (d + shackle[2].s) / fd0);
          s_cat = (d + shackle[2].s) * (1 + 2 * (fd0 / (spm_chain_sub * g * (d + shackle[2].s)))) ^ 0.5;
          x_lnk_cat[1] = x_cat;
          x_lnk_plot_o = L_chain - s_cat + x_lnk_cat[1];
          x_lnk_plot[1] = x_lnk_plot_o - (x_lnk_plot_o - shackle[1].s);
          z_lnk_cat[1] = shackle[2].s + d + a;
          z_lnk_plot[1] = (-(d + a)) + z_lnk_cat[1];
          for i in 2:n_lnk loop
            if x_lnk_cat[i - 1] > 0 then
              x_lnk_cat[i] = a * Modelica.Math.asinh((s_cat - (i - 1) * lnk_l) / a);
              if x_lnk_cat[i] > 0 then
                x_lnk_plot[i] = L_chain - s_cat + x_lnk_cat[i] - (x_lnk_plot_o - shackle[1].s);
                z_lnk_cat[i] = a * cosh(x_lnk_cat[i] / a);
                z_lnk_plot[i] = (-(d + a)) + z_lnk_cat[i];
              else
                x_lnk_plot[i] = x_lnk_plot[i - 1] - sqrt(lnk_l ^ 2 - (d + z_lnk_plot[i - 1]) ^ 2);
                z_lnk_cat[i] = 0;
                z_lnk_plot[i] = -d;
              end if;
            else
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            end if;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          shackle[1].f = fd0;
          shackle[2].f = s_cat * spm_chain_sub * g;
        else
          fd0 = 0;
          a = 0;
          z_cat = 0;
          psi_cat = 0;
          x_cat = 0;
          s_cat = 0;
          x_lnk_cat[1] = 0;
          x_lnk_plot_o = 0;
          x_lnk_plot[1] = n_lnk * lnk_l;
          z_lnk_cat[1] = 0;
          z_lnk_plot[1] = -d;
          for i in 2:n_lnk loop
            x_lnk_cat[i] = 0;
            x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
            z_lnk_cat[i] = 0;
            z_lnk_plot[i] = -d;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          shackle[1].f = 0;
          shackle[2].f = 0;
        end if;
        emdc.xinit = X0 - shackle[2].s;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -92}, fillColor = {170, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 8}, {100, -8}}), Rectangle(origin = {0, -18}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 68}, {100, -68}}), Rectangle(origin = {-87, -78}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-9, 8}, {9, -8}}), Line(origin = {1.26, -16.13}, points = {{-78.9938, -64.1437}, {-68.9938, -70.1437}, {-8.9938, -70.1437}, {33.0062, -36.1437}, {65.0062, 15.8563}, {79.0062, 63.8563}}, color = {255, 0, 0}, thickness = 0.5, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-06, Interval = 0.5));
      end CatenaryMooring_Mf0;

      model CatenaryMooring_MfC
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_b shackle[2] annotation(
          Placement(transformation(extent = {{-20, 80}, {20, 120}})));
        OceanEngineering.Connectors.EnvironmentMooringDataConnector emdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{80, -20}, {120, 20}}, rotation = 0), iconTransformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        parameter Integer nOmega = 100 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        constant Real g = Modelica.Constants.g_n;
        constant Real pi = Modelica.Constants.pi;
        parameter Real rho_w = 1025;
        parameter Real d = 100;
        parameter Real lnk_l = 3;
        parameter Integer n_lnk = 50;
        parameter Real L_chain = lnk_l * n_lnk;
        parameter Real d_chain = 0.022;
        parameter Real rho_mat_chain = 7800;
        parameter Real spm_chain = 10;
        parameter Real spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Real Xmax = sqrt(L_chain ^ 2 - d ^ 2);
        parameter Real X0 = L_chain - d;
        parameter Real xmax = Xmax - 10;
        parameter Real x[:] = 0:0.1:xmax;
        parameter Real Th[:] = OceanEngineering.Functions.catThIterator(L_chain, x, d, spm_chain_sub);
        parameter Real X[:] = OceanEngineering.Functions.catXIterator(L_chain, d, spm_chain_sub, Th, x);
        parameter Real Cdn = 1;
        parameter Real Cdt = 0;
        parameter Real Cman = 1;
        parameter Real Cmat = 0.0;
        parameter Real Cmn = 1 + Cman;
        parameter Real Cmt = 1 + Cmat;
        parameter Real D = 1.8 * d_chain;
        parameter Real Ai = rho_w * pi * D ^ 2 / 4;
        parameter Real Ad = 0.5 * rho_w * D;
        Real x_lnk_cat[n_lnk + 1];
        Real x_lnk_plot_o;
        Real z_lnk_cat[n_lnk + 1];
        Real x_lnk_plot[n_lnk + 1];
        Real z_lnk_plot[n_lnk + 1];
        Real xl[n_lnk];
        Real zl[n_lnk];
        Real vxl[n_lnk];
        Real vzl[n_lnk];
        Real psil[n_lnk];
        Real vl[n_lnk];
        Real vl_ang[n_lnk];
        Real vln[n_lnk];
        Real vlt[n_lnk];
        Real fd0;
        Real a;
        Real z_cat;
        Real x_cat;
        Real s_cat;
        Real SSE_Xl[n_lnk];
        Real Uc[n_lnk];
        Real Ucn[n_lnk];
        Real Uct[n_lnk];
        Real Mfni[n_lnk];
        Real Mfti[n_lnk];
        Real Mfxi[n_lnk];
        Real Mfzi[n_lnk];
        Real Mfx;
        Real Mfz;
        parameter Real Trmp = 50;
        Real Mfxa;
        Real Mfza;
      equation
        if shackle[1].s <= X0 - shackle[2].s then
          fd0 = 0;
          a = 0;
          z_cat = d + shackle[2].s;
          x_cat = X0 - shackle[2].s;
          s_cat = d + shackle[2].s;
          x_lnk_cat[1] = 0;
          x_lnk_plot_o = shackle[1].s;
          x_lnk_plot[1] = shackle[1].s;
          z_lnk_cat[1] = 0;
          z_lnk_plot[1] = shackle[2].s;
          for i in 2:n_lnk loop
            if z_lnk_plot[i - 1] > (-(d - lnk_l)) then
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = shackle[1].s;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = z_lnk_plot[i - 1] - lnk_l;
            elseif z_lnk_plot[i - 1] > (-d) then
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - sqrt(lnk_l ^ 2 - (d + z_lnk_plot[i - 1]) ^ 2);
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            else
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            end if;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          for i in 1:n_lnk loop
            xl[i] = (x_lnk_plot[i] + x_lnk_plot[i + 1]) / 2;
            zl[i] = (z_lnk_plot[i] + z_lnk_plot[i + 1]) / 2;
          end for;
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 then
              vxl[i] = if noEvent(der(xl[i])) then der(xl[i]) else 0;
              vzl[i] = if noEvent(der(zl[i])) then der(zl[i]) else 0;
            else
              vxl[i] = 0;
              vzl[i] = 0;
            end if;
          end for;
// Mf Calculation loop 1
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 * lnk_l then
              psil[i] = if noEvent(abs(x_lnk_plot[i] - x_lnk_plot[i + 1]) > 0) then atan((z_lnk_plot[i] - z_lnk_plot[i + 1]) / (x_lnk_plot[i] - x_lnk_plot[i + 1])) else 0;
              SSE_Xl[i] = OceanEngineering.Functions.sSE_X_cat(time, emdc.omega, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i]);
              vl[i] = sqrt(vxl[i] ^ 2 + vzl[i] ^ 2);
              vl_ang[i] = if noEvent(abs(vxl[i]) > 0) then atan(vzl[i] / vxl[i]) else pi / 2;
              vln[i] = -vl[i] * sin(psil[i] - vl_ang[i]);
              vlt[i] = vl[i] * cos(psil[i] - vl_ang[i]);
              Uc[i] = OceanEngineering.Functions.linearInterpolatorSV(emdc.zcg, emdc.Ucg, zl[i] + SSE_Xl[i]);
              Ucn[i] = Uc[i] * sin(psil[i]);
              Uct[i] = Uc[i] * cos(psil[i]);
              Mfni[i] = Cdn * Ad * abs(Ucn[i] - vln[i]) * (Ucn[i] - vln[i]) * lnk_l;
              Mfti[i] = Cdt * Ad * abs(Uct[i] - vlt[i]) * (Uct[i] - vlt[i]) * lnk_l;
              Mfxi[i] = Mfni[i] * sin(psil[i]) + Mfti[i] * cos(psil[i]);
              Mfzi[i] = Mfni[i] * cos(psil[i]) + Mfti[i] * sin(psil[i]);
            else
              psil[i] = 0;
              SSE_Xl[i] = 0;
              vl[i] = 0;
              vl_ang[i] = 0;
              vln[i] = 0;
              vlt[i] = 0;
              Uc[i] = 0;
              Ucn[i] = 0;
              Uct[i] = 0;
              Mfni[i] = 0;
              Mfti[i] = 0;
              Mfxi[i] = 0;
              Mfzi[i] = 0;
            end if;
          end for;
          Mfx = sum(Mfxi);
          Mfz = sum(Mfzi);
          if time < Trmp then
            Mfxa = time / Trmp * Mfx;
            Mfza = time / Trmp * Mfz;
          else
            Mfxa = Mfx;
            Mfza = Mfz;
          end if;
          shackle[1].f = fd0 - Mfxa - 10;
          shackle[2].f = s_cat * spm_chain_sub * g + Mfza;
//MAIN IF cond 2
        elseif shackle[1].s < sqrt(L_chain ^ 2 - d ^ 2) then
          fd0 = OceanEngineering.Functions.linearInterpolatorSV(X, Th, shackle[1].s);
          a = fd0 / (spm_chain_sub * g);
          z_cat = a + d + shackle[2].s;
          x_cat = fd0 / (spm_chain_sub * g) * Modelica.Math.acosh(1 + spm_chain_sub * g * (d + shackle[2].s) / fd0);
          s_cat = (d + shackle[2].s) * (1 + 2 * (fd0 / (spm_chain_sub * g * (d + shackle[2].s)))) ^ 0.5;
          x_lnk_cat[1] = x_cat;
          x_lnk_plot_o = L_chain - s_cat + x_lnk_cat[1];
          x_lnk_plot[1] = x_lnk_plot_o - (x_lnk_plot_o - shackle[1].s);
          z_lnk_cat[1] = shackle[2].s + d + a;
          z_lnk_plot[1] = (-(d + a)) + z_lnk_cat[1];
          for i in 2:n_lnk loop
            if x_lnk_cat[i - 1] > 0 then
              x_lnk_cat[i] = a * Modelica.Math.asinh((s_cat - (i - 1) * lnk_l) / a);
              if x_lnk_cat[i] > 0 then
                x_lnk_plot[i] = L_chain - s_cat + x_lnk_cat[i] - (x_lnk_plot_o - shackle[1].s);
                z_lnk_cat[i] = a * cosh(x_lnk_cat[i] / a);
                z_lnk_plot[i] = (-(d + a)) + z_lnk_cat[i];
              else
                x_lnk_plot[i] = x_lnk_plot[i - 1] - sqrt(lnk_l ^ 2 - (d + z_lnk_plot[i - 1]) ^ 2);
                z_lnk_cat[i] = 0;
                z_lnk_plot[i] = -d;
              end if;
            else
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            end if;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          for i in 1:n_lnk loop
            xl[i] = (x_lnk_plot[i] + x_lnk_plot[i + 1]) / 2;
            zl[i] = (z_lnk_plot[i] + z_lnk_plot[i + 1]) / 2;
          end for;
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 then
              vxl[i] = if noEvent(der(xl[i])) then der(xl[i]) else 0;
              vzl[i] = if noEvent(der(zl[i])) then der(zl[i]) else 0;
            else
              vxl[i] = 0;
              vzl[i] = 0;
            end if;
          end for;
//Mf calculation loop 2
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 * lnk_l then
              psil[i] = if noEvent(abs(x_lnk_plot[i] - x_lnk_plot[i + 1]) > 0) then atan((z_lnk_plot[i] - z_lnk_plot[i + 1]) / (x_lnk_plot[i] - x_lnk_plot[i + 1])) else 0;
              SSE_Xl[i] = OceanEngineering.Functions.sSE_X_cat(time, emdc.omega, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i]);
              vl[i] = sqrt(vxl[i] ^ 2 + vzl[i] ^ 2);
              vl_ang[i] = if noEvent(abs(vxl[i]) > 0) then atan(vzl[i] / vxl[i]) else pi / 2;
              vln[i] = -vl[i] * sin(psil[i] - vl_ang[i]);
              vlt[i] = vl[i] * cos(psil[i] - vl_ang[i]);
              Uc[i] = OceanEngineering.Functions.linearInterpolatorSV(emdc.zcg, emdc.Ucg, zl[i] + SSE_Xl[i]);
              Ucn[i] = Uc[i] * sin(psil[i]);
              Uct[i] = Uc[i] * cos(psil[i]);
              Mfni[i] = Cdn * Ad * abs(Ucn[i] - vln[i]) * (Ucn[i] - vln[i]) * lnk_l;
              Mfti[i] = Cdt * Ad * abs(Uct[i] - vlt[i]) * (Uct[i] - vlt[i]) * lnk_l;
              Mfxi[i] = Mfni[i] * sin(psil[i]) + Mfti[i] * cos(psil[i]);
              Mfzi[i] = Mfni[i] * cos(psil[i]) + Mfti[i] * sin(psil[i]);
            else
              psil[i] = 0;
              SSE_Xl[i] = 0;
              vl[i] = 0;
              vl_ang[i] = 0;
              vln[i] = 0;
              vlt[i] = 0;
              Uc[i] = 0;
              Ucn[i] = 0;
              Uct[i] = 0;
              Mfni[i] = 0;
              Mfti[i] = 0;
              Mfxi[i] = 0;
              Mfzi[i] = 0;
            end if;
          end for;
          Mfx = sum(Mfxi);
          Mfz = sum(Mfzi);
          if time < Trmp then
            Mfxa = time / Trmp * Mfx;
            Mfza = time / Trmp * Mfz;
          else
            Mfxa = Mfx;
            Mfza = Mfz;
          end if;
          shackle[1].f = fd0 - Mfxa;
          shackle[2].f = s_cat * spm_chain_sub * g + Mfza;
//MAIN IF cond 3
        else
          fd0 = 0;
          a = 0;
          z_cat = 0;
          x_cat = 0;
          s_cat = 0;
          x_lnk_cat[1] = 0;
          x_lnk_plot_o = 0;
          x_lnk_plot[1] = lnk_l * n_lnk;
          z_lnk_cat[1] = 0;
          z_lnk_plot[1] = -d;
          for i in 2:n_lnk loop
            x_lnk_cat[i] = 0;
            x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
            z_lnk_cat[i] = 0;
            z_lnk_plot[i] = -d;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          for i in 1:n_lnk loop
            xl[i] = 0;
            zl[i] = 0;
          end for;
          for i in 1:n_lnk loop
            vxl[i] = 0;
            vzl[i] = 0;
          end for;
//Mf calculation loop 2
          for i in 1:n_lnk loop
            psil[i] = 0;
            SSE_Xl[i] = 0;
            vl[i] = 0;
            vl_ang[i] = 0;
            vln[i] = 0;
            vlt[i] = 0;
            Uc[i] = 0;
            Ucn[i] = 0;
            Uct[i] = 0;
            Mfni[i] = 0;
            Mfti[i] = 0;
            Mfxi[i] = 0;
            Mfzi[i] = 0;
          end for;
          Mfx = 0;
          Mfz = 0;
          Mfxa = 0;
          Mfza = 0;
          shackle[1].f = 0;
          shackle[2].f = 0;
        end if;
        emdc.xinit = X0 - shackle[2].s;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -92}, fillColor = {170, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 8}, {100, -8}}), Rectangle(origin = {0, -18}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 68}, {100, -68}}), Rectangle(origin = {-87, -78}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-9, 8}, {9, -8}}), Line(origin = {0.99, -15.86}, points = {{-78.9938, -64.1437}, {-68.9938, -70.1437}, {-8.9938, -70.1437}, {33.0062, -36.1437}, {65.0062, 15.8563}, {79.0062, 63.8563}}, color = {255, 0, 0}, thickness = 0.5, smooth = Smooth.Bezier), Line(origin = {-22.3805, 26.2699}, points = {{-34, 0}, {34, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-21.2879, 6.5784}, points = {{-34, 0}, {20, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-20.7352, -10.6838}, points = {{-34, 0}, {4, 0}}, arrow = {Arrow.None, Arrow.Filled}), Text(origin = {-33, -27}, extent = {{-21, 7}, {21, -7}}, textString = "Current")}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-6, Interval = 0.5));
      end CatenaryMooring_MfC;


      model CatenaryMooring_MfCW
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_b shackle[2] annotation(
          Placement(transformation(extent = {{-20, 80}, {20, 120}})));
        OceanEngineering.Connectors.EnvironmentMooringDataConnector emdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{80, -20}, {120, 20}}, rotation = 0), iconTransformation(extent = {{-120, -20}, {-80, 20}}, rotation = 0)));
        parameter Integer nOmega = 100 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        constant Real g = Modelica.Constants.g_n;
        constant Real pi = Modelica.Constants.pi;
        parameter Real rho_w = 1025;
        parameter Real d = 100;
        parameter Real lnk_l = 3;
        parameter Integer n_lnk = 50;
        parameter Real L_chain = lnk_l * n_lnk;
        parameter Real d_chain = 0.022;
        parameter Real rho_mat_chain = 7800;
        parameter Real spm_chain = 10;
        parameter Real Cdn = 1;
        parameter Real Cdt = 0;
        parameter Real Cman = 1;
        parameter Real Cmat = 0;
        parameter Real Trmp = 50;
        parameter Real spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Real Xmax = sqrt(L_chain ^ 2 - d ^ 2);
        parameter Real X0 = L_chain - d;
        parameter Real xmax = Xmax - 10;
        parameter Real x[:] = 0:0.1:xmax;
        parameter Real Th[:] = OceanEngineering.Functions.catThIterator(L_chain, x, d, spm_chain_sub);
        parameter Real X[:] = OceanEngineering.Functions.catXIterator(L_chain, d, spm_chain_sub, Th, x);
        parameter Real Cmn = 1 + Cman;
        parameter Real Cmt = 1 + Cmat;
        parameter Real D = 1.8 * d_chain;
        parameter Real Ai = rho_w * pi * D ^ 2 / 4;
        parameter Real Ad = 0.5 * rho_w * D;
        Real x_lnk_cat[n_lnk + 1];
        Real x_lnk_plot_o;
        Real z_lnk_cat[n_lnk + 1];
        Real x_lnk_plot[n_lnk + 1];
        Real z_lnk_plot[n_lnk + 1];
        Real xl[n_lnk];
        Real zl[n_lnk];
        Real vxl[n_lnk];
        Real axl[n_lnk];
        Real vzl[n_lnk];
        Real azl[n_lnk];
        Real psil[n_lnk];
        Real vl[n_lnk];
        Real vl_ang[n_lnk];
        Real vln[n_lnk];
        Real vlt[n_lnk];
        Real al[n_lnk];
        Real al_ang[n_lnk];
        Real aln[n_lnk];
        Real alt[n_lnk];
        Real ul[n_lnk];
        Real wl[n_lnk];
        Real vw[n_lnk];
        Real vw_ang[n_lnk];
        Real vwn[n_lnk];
        Real vwt[n_lnk];
        Real aul[n_lnk];
        Real awl[n_lnk];
        Real aw[n_lnk];
        Real aw_ang[n_lnk];
        Real awn[n_lnk];
        Real awt[n_lnk];
        Real MCman[n_lnk];
        Real MCmat[n_lnk];
        Real fd0;
        Real a;
        Real z_cat;
        Real x_cat;
        Real s_cat;
        Real SSE_Xl[n_lnk];
        Real Uc[n_lnk];
        Real Ucn[n_lnk];
        Real Uct[n_lnk];
        Real Mfni[n_lnk];
        Real Mfti[n_lnk];
        Real Mfxi[n_lnk];
        Real Mfzi[n_lnk];
        Real Mfx;
        Real Mfz;
        Real Mfxa;
        Real Mfza;
      equation
        if shackle[1].s <= X0 - shackle[2].s then
          fd0 = 0;
          a = 0;
          z_cat = d + shackle[2].s;
          x_cat = X0 - shackle[2].s;
          s_cat = d + shackle[2].s;
          x_lnk_cat[1] = 0;
          x_lnk_plot_o = shackle[1].s;
          x_lnk_plot[1] = shackle[1].s;
          z_lnk_cat[1] = 0;
          z_lnk_plot[1] = shackle[2].s;
          for i in 2:n_lnk loop
            if z_lnk_plot[i - 1] > (-(d - lnk_l)) then
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = shackle[1].s;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = z_lnk_plot[i - 1] - lnk_l;
            elseif z_lnk_plot[i - 1] > (-d) then
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - sqrt(lnk_l ^ 2 - (d + z_lnk_plot[i - 1]) ^ 2);
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            else
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            end if;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          for i in 1:n_lnk loop
            xl[i] = (x_lnk_plot[i] + x_lnk_plot[i + 1]) / 2;
            zl[i] = (z_lnk_plot[i] + z_lnk_plot[i + 1]) / 2;
          end for;
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 then
              vxl[i] = if noEvent(der(xl[i])) then der(xl[i]) else 0;
              axl[i] = if noEvent(der(vxl[i])) then der(vxl[i]) else 0;
              vzl[i] = if noEvent(der(zl[i])) then der(zl[i]) else 0;
              azl[i] = if noEvent(der(vzl[i])) then der(vzl[i]) else 0;
            else
              vxl[i] = 0;
              axl[i] = 0;
              vzl[i] = 0;
              azl[i] = 0;
            end if;
          end for;
// Mf Calculation loop 1
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 * lnk_l then
              psil[i] = if noEvent(abs(x_lnk_plot[i] - x_lnk_plot[i + 1]) > 0) then atan((z_lnk_plot[i] - z_lnk_plot[i + 1]) / (x_lnk_plot[i] - x_lnk_plot[i + 1])) else 0;
              SSE_Xl[i] = OceanEngineering.Functions.sSE_X_cat(time, emdc.omega, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i]);
              vl[i] = sqrt(vxl[i] ^ 2 + vzl[i] ^ 2);
              vl_ang[i] = if noEvent(abs(vxl[i]) > 0) then atan(vzl[i] / vxl[i]) else pi / 2;
              vln[i] = -vl[i] * sin(psil[i] - vl_ang[i]);
              vlt[i] = vl[i] * cos(psil[i] - vl_ang[i]);
              al[i] = sqrt(axl[i] ^ 2 + azl[i] ^ 2);
              al_ang[i] = if noEvent(abs(axl[i]) > 0) then atan(azl[i] / axl[i]) else pi / 2;
              aln[i] = -al[i] * sin(psil[i] - al_ang[i]);
              alt[i] = al[i] * cos(psil[i] - al_ang[i]);
              ul[i] = OceanEngineering.Functions.wave_uCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              wl[i] = OceanEngineering.Functions.wave_wCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              vw[i] = sqrt(ul[i] ^ 2 + wl[i] ^ 2);
              vw_ang[i] = if noEvent(abs(ul[i]) > 0) then atan(wl[i] / ul[i]) else pi / 2;
              vwn[i] = -vw[i] * sin(psil[i] - vw_ang[i]);
              vwt[i] = vw[i] * cos(psil[i] - vw_ang[i]);
              aul[i] = OceanEngineering.Functions.wave_auCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              awl[i] = OceanEngineering.Functions.wave_awCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              aw[i] = sqrt(aul[i] ^ 2 + awl[i] ^ 2);
              aw_ang[i] = if noEvent(abs(aul[i]) > 0) then atan(awl[i] / aul[i]) else pi / 2;
              awn[i] = -aw[i] * sin(psil[i] - aw_ang[i]);
              awt[i] = aw[i] * cos(psil[i] - aw_ang[i]);
              if abs(vw[i]) > 0 or abs(aw[i]) > 0 then
                MCman[i] = Cman * Ai * aln[i];
                MCmat[i] = Cmat * Ai * alt[i];
              else
                MCman[i] = 0;
                MCmat[i] = 0;
              end if;
              Uc[i] = OceanEngineering.Functions.linearInterpolatorSV(emdc.zcg, emdc.Ucg, zl[i] + SSE_Xl[i]);
              Ucn[i] = Uc[i] * sin(psil[i]);
              Uct[i] = Uc[i] * cos(psil[i]);
              if zl[i] > SSE_Xl[i] then
                Mfni[i] = 0;
                Mfti[i] = 0;
                Mfxi[i] = 0;
                Mfzi[i] = 0;
              else
                Mfni[i] = (Cmn * Ai * awn[i] - MCman[i] + Cdn * Ad * abs(vwn[i] + Ucn[i] - vln[i]) * (vwn[i] + Ucn[i] - vln[i])) * lnk_l;
                Mfti[i] = (Cmt * Ai * awt[i] - MCmat[i] + Cdt * Ad * abs(vwt[i] + Uct[i] - vlt[i]) * (vwt[i] + Uct[i] - vlt[i])) * lnk_l;
                Mfxi[i] = Mfni[i] * sin(psil[i]) + Mfti[i] * cos(psil[i]);
                Mfzi[i] = Mfni[i] * cos(psil[i]) + Mfti[i] * sin(psil[i]);
              end if;
            else
              psil[i] = 0;
              SSE_Xl[i] = 0;
              vl[i] = 0;
              vl_ang[i] = 0;
              vln[i] = 0;
              vlt[i] = 0;
              al[i] = 0;
              al_ang[i] = 0;
              aln[i] = 0;
              alt[i] = 0;
              ul[i] = 0;
              wl[i] = 0;
              vw[i] = 0;
              vw_ang[i] = 0;
              vwn[i] = 0;
              vwt[i] = 0;
              aul[i] = 0;
              awl[i] = 0;
              aw[i] = 0;
              aw_ang[i] = 0;
              awn[i] = 0;
              awt[i] = 0;
              MCman[i] = 0;
              MCmat[i] = 0;
              Uc[i] = 0;
              Ucn[i] = 0;
              Uct[i] = 0;
              Mfni[i] = 0;
              Mfti[i] = 0;
              Mfxi[i] = 0;
              Mfzi[i] = 0;
            end if;
          end for;
          Mfx = sum(Mfxi);
          Mfz = sum(Mfzi);
          if time < Trmp then
            Mfxa = time / Trmp * Mfx;
            Mfza = time / Trmp * Mfz;
          else
            Mfxa = Mfx;
            Mfza = Mfz;
          end if;
          shackle[1].f = fd0 - Mfxa - 10;
          shackle[2].f = s_cat * spm_chain_sub * g + Mfza;
//MAIN IF Cond2
        elseif shackle[1].s < sqrt(L_chain ^ 2 - d ^ 2) then
          fd0 = OceanEngineering.Functions.linearInterpolatorSV(X, Th, shackle[1].s);
          a = fd0 / (spm_chain_sub * g);
          z_cat = a + d + shackle[2].s;
          x_cat = fd0 / (spm_chain_sub * g) * Modelica.Math.acosh(1 + spm_chain_sub * g * (d + shackle[2].s) / fd0);
          s_cat = (d + shackle[2].s) * (1 + 2 * (fd0 / (spm_chain_sub * g * (d + shackle[2].s)))) ^ 0.5;
          x_lnk_cat[1] = x_cat;
          x_lnk_plot_o = L_chain - s_cat + x_lnk_cat[1];
          x_lnk_plot[1] = x_lnk_plot_o - (x_lnk_plot_o - shackle[1].s);
          z_lnk_cat[1] = shackle[2].s + d + a;
          z_lnk_plot[1] = (-(d + a)) + z_lnk_cat[1];
          for i in 2:n_lnk loop
            if x_lnk_cat[i - 1] > 0 then
              x_lnk_cat[i] = a * Modelica.Math.asinh((s_cat - (i - 1) * lnk_l) / a);
              if x_lnk_cat[i] > 0 then
                x_lnk_plot[i] = L_chain - s_cat + x_lnk_cat[i] - (x_lnk_plot_o - shackle[1].s);
                z_lnk_cat[i] = a * cosh(x_lnk_cat[i] / a);
                z_lnk_plot[i] = (-(d + a)) + z_lnk_cat[i];
              else
                x_lnk_plot[i] = x_lnk_plot[i - 1] - sqrt(lnk_l ^ 2 - (d + z_lnk_plot[i - 1]) ^ 2);
                z_lnk_cat[i] = 0;
                z_lnk_plot[i] = -d;
              end if;
            else
              x_lnk_cat[i] = 0;
              x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
              z_lnk_cat[i] = 0;
              z_lnk_plot[i] = -d;
            end if;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          for i in 1:n_lnk loop
            xl[i] = (x_lnk_plot[i] + x_lnk_plot[i + 1]) / 2;
            zl[i] = (z_lnk_plot[i] + z_lnk_plot[i + 1]) / 2;
          end for;
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 then
              vxl[i] = if noEvent(der(xl[i])) then der(xl[i]) else 0;
              axl[i] = if noEvent(der(vxl[i])) then der(vxl[i]) else 0;
              vzl[i] = if noEvent(der(zl[i])) then der(zl[i]) else 0;
              azl[i] = if noEvent(der(vzl[i])) then der(vzl[i]) else 0;
            else
              vxl[i] = 0;
              axl[i] = 0;
              vzl[i] = 0;
              azl[i] = 0;
            end if;
          end for;
//Mf calculation loop 2
          for i in 1:n_lnk loop
            if zl[i] > (-d) + 1 * lnk_l then
              psil[i] = if noEvent(abs(x_lnk_plot[i] - x_lnk_plot[i + 1]) > 0) then atan((z_lnk_plot[i] - z_lnk_plot[i + 1]) / (x_lnk_plot[i] - x_lnk_plot[i + 1])) else 0;
              SSE_Xl[i] = OceanEngineering.Functions.sSE_X_cat(time, emdc.omega, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i]);
              vl[i] = sqrt(vxl[i] ^ 2 + vzl[i] ^ 2);
              vl_ang[i] = if noEvent(abs(vxl[i]) > 0) then atan(vzl[i] / vxl[i]) else pi / 2;
              vln[i] = -vl[i] * sin(psil[i] - vl_ang[i]);
              vlt[i] = vl[i] * cos(psil[i] - vl_ang[i]);
              al[i] = sqrt(axl[i] ^ 2 + azl[i] ^ 2);
              al_ang[i] = if noEvent(abs(axl[i]) > 0) then atan(azl[i] / axl[i]) else pi / 2;
              aln[i] = -al[i] * sin(psil[i] - al_ang[i]);
              alt[i] = al[i] * cos(psil[i] - al_ang[i]);
              ul[i] = OceanEngineering.Functions.wave_uCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              wl[i] = OceanEngineering.Functions.wave_wCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              vw[i] = sqrt(ul[i] ^ 2 + wl[i] ^ 2);
              vw_ang[i] = if noEvent(abs(ul[i]) > 0) then atan(wl[i] / ul[i]) else pi / 2;
              vwn[i] = -vw[i] * sin(psil[i] - vw_ang[i]);
              vwt[i] = vw[i] * cos(psil[i] - vw_ang[i]);
              aul[i] = OceanEngineering.Functions.wave_auCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              awl[i] = OceanEngineering.Functions.wave_awCalculator(time, d, emdc.omega, emdc.T, emdc.k, emdc.zeta0i, emdc.epsilon, xl[i], zl[i], SSE_Xl[i]);
              aw[i] = sqrt(aul[i] ^ 2 + awl[i] ^ 2);
              aw_ang[i] = if noEvent(abs(aul[i]) > 0) then atan(awl[i] / aul[i]) else pi / 2;
              awn[i] = -aw[i] * sin(psil[i] - aw_ang[i]);
              awt[i] = aw[i] * cos(psil[i] - aw_ang[i]);
              if abs(vw[i]) > 0 or abs(aw[i]) > 0 then
                MCman[i] = Cman * Ai * aln[i];
                MCmat[i] = Cmat * Ai * alt[i];
              else
                MCman[i] = 0;
                MCmat[i] = 0;
              end if;
              Uc[i] = OceanEngineering.Functions.linearInterpolatorSV(emdc.zcg, emdc.Ucg, zl[i] + SSE_Xl[i]);
              Ucn[i] = Uc[i] * sin(psil[i]);
              Uct[i] = Uc[i] * cos(psil[i]);
              if zl[i] > SSE_Xl[i] then
                Mfni[i] = 0;
                Mfti[i] = 0;
                Mfxi[i] = 0;
                Mfzi[i] = 0;
              else
                Mfni[i] = (Cmn * Ai * awn[i] - MCman[i] + Cdn * Ad * abs(vwn[i] + Ucn[i] - vln[i]) * (vwn[i] + Ucn[i] - vln[i])) * lnk_l;
                Mfti[i] = (Cmt * Ai * awt[i] - MCmat[i] + Cdt * Ad * abs(vwt[i] + Uct[i] - vlt[i]) * (vwt[i] + Uct[i] - vlt[i])) * lnk_l;
                Mfxi[i] = Mfni[i] * sin(psil[i]) + Mfti[i] * cos(psil[i]);
                Mfzi[i] = Mfni[i] * cos(psil[i]) + Mfti[i] * sin(psil[i]);
              end if;
            else
              psil[i] = 0;
              SSE_Xl[i] = 0;
              vl[i] = 0;
              vl_ang[i] = 0;
              vln[i] = 0;
              vlt[i] = 0;
              al[i] = 0;
              al_ang[i] = 0;
              aln[i] = 0;
              alt[i] = 0;
              ul[i] = 0;
              wl[i] = 0;
              vw[i] = 0;
              vw_ang[i] = 0;
              vwn[i] = 0;
              vwt[i] = 0;
              aul[i] = 0;
              awl[i] = 0;
              aw[i] = 0;
              aw_ang[i] = 0;
              awn[i] = 0;
              awt[i] = 0;
              MCman[i] = 0;
              MCmat[i] = 0;
              Uc[i] = 0;
              Ucn[i] = 0;
              Uct[i] = 0;
              Mfni[i] = 0;
              Mfti[i] = 0;
              Mfxi[i] = 0;
              Mfzi[i] = 0;
            end if;
          end for;
          Mfx = sum(Mfxi);
          Mfz = sum(Mfzi);
          if time < Trmp then
            Mfxa = time / Trmp * Mfx;
            Mfza = time / Trmp * Mfz;
          else
            Mfxa = Mfx;
            Mfza = Mfz;
          end if;
          shackle[1].f = fd0 - Mfxa;
          shackle[2].f = s_cat * spm_chain_sub * g + Mfza;
//MAIN IF Cond3
        else
          fd0 = 0;
          a = 0;
          z_cat = 0;
          x_cat = 0;
          s_cat = 0;
          x_lnk_cat[1] = 0;
          x_lnk_plot_o = 0;
          x_lnk_plot[1] = n_lnk * lnk_l;
          z_lnk_cat[1] = 0;
          z_lnk_plot[1] = -d;
          for i in 2:n_lnk loop
            x_lnk_cat[i] = 0;
            x_lnk_plot[i] = x_lnk_plot[i - 1] - lnk_l;
            z_lnk_cat[i] = 0;
            z_lnk_plot[i] = -d;
          end for;
          x_lnk_cat[n_lnk + 1] = 0;
          x_lnk_plot[n_lnk + 1] = 0;
          z_lnk_cat[n_lnk + 1] = 0;
          z_lnk_plot[n_lnk + 1] = -d;
          for i in 1:n_lnk loop
            xl[i] = 0;
            zl[i] = 0;
          end for;
          for i in 1:n_lnk loop
            vxl[i] = 0;
            axl[i] = 0;
            vzl[i] = 0;
            azl[i] = 0;
          end for;
//Mf calculation loop 2
          for i in 1:n_lnk loop
            psil[i] = 0;
            SSE_Xl[i] = 0;
            vl[i] = 0;
            vl_ang[i] = 0;
            vln[i] = 0;
            vlt[i] = 0;
            al[i] = 0;
            al_ang[i] = 0;
            aln[i] = 0;
            alt[i] = 0;
            ul[i] = 0;
            wl[i] = 0;
            vw[i] = 0;
            vw_ang[i] = 0;
            vwn[i] = 0;
            vwt[i] = 0;
            aul[i] = 0;
            awl[i] = 0;
            aw[i] = 0;
            aw_ang[i] = 0;
            awn[i] = 0;
            awt[i] = 0;
            MCman[i] = 0;
            MCmat[i] = 0;
            Uc[i] = 0;
            Ucn[i] = 0;
            Uct[i] = 0;
            Mfni[i] = 0;
            Mfti[i] = 0;
            Mfxi[i] = 0;
            Mfzi[i] = 0;
          end for;
          Mfx = 0;
          Mfz = 0;
          Mfxa = Mfx;
          Mfza = Mfz;
          shackle[1].f = 0;
          shackle[2].f = 0;
        end if;
        emdc.xinit = X0 - shackle[2].s;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -92}, fillColor = {170, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 8}, {100, -8}}), Rectangle(origin = {0, -18}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 68}, {100, -68}}), Rectangle(origin = {-87, -78}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-9, 8}, {9, -8}}), Line(origin = {0.99, -15.86}, points = {{-78.9938, -64.1437}, {-68.9938, -70.1437}, {-8.9938, -70.1437}, {33.0062, -36.1437}, {65.0062, 15.8563}, {79.0062, 63.8563}}, color = {255, 0, 0}, thickness = 1, smooth = Smooth.Bezier), Line(origin = {-30, 22}, points = {{-32, 0}, {32, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-29.4473, 5.00771}, points = {{-32, 0}, {20, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-28.8946, -11.7147}, points = {{-32, 0}, {8, 0}}, arrow = {Arrow.None, Arrow.Filled}), Text(origin = {-44, -28}, extent = {{-18, 8}, {18, -8}}, textString = "Current"), Text(origin = {-70, 70}, extent = {{-18, 8}, {18, -8}}, textString = "Waves"), Polygon(origin = {-15, 49}, lineColor = {0, 170, 255}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0, points = {{-97, -7}, {-69, 11}, {-47, -3}, {-21, 13}, {7, -3}, {35, 17}, {53, -3}, {101, 11}, {127, -7}, {-97, -7}}, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.5));
      end CatenaryMooring_MfCW;


    end Moorings;
  end Components;

  package MiscellaneousComponents
    model DampedSpring_Vertical
      parameter StateSelect stateSelect = StateSelect.prefer "Priority to use s_rel and v_rel as states" annotation(
        HideResult = true,
        Dialog(tab = "Advanced"));
      parameter Modelica.SIunits.Distance s_nominal = 1e-4 "Nominal value of s_rel (used for scaling)" annotation(
        Dialog(tab = "Advanced"));
      Modelica.SIunits.Position s_rel(start = 0.9019, stateSelect = stateSelect, nominal = s_nominal) "Relative distance (= flange_b.s - flange_a.s)";
      Modelica.SIunits.Velocity v_rel(start = 0, stateSelect = stateSelect) "Relative velocity (= der(s_rel))";
      Modelica.SIunits.Force f "Forces between flanges (= flange_b.f)";
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a "Left flange of compliant 1-dim. translational component" annotation(
        Placement(transformation(extent = {{-10, -110}, {10, -90}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Right flange of compliant 1-dim. translational component" annotation(
        Placement(transformation(extent = {{-10, 90}, {10, 110}})));
      parameter Modelica.SIunits.TranslationalSpringConstant c(final min = 0, start = 5000) "Spring constant";
      parameter Modelica.SIunits.TranslationalDampingConstant d(final min = 0, start = 50) "Damping constant";
      parameter Modelica.SIunits.Position s_rel0 = 0.9019 "Unstretched spring length";
      extends Modelica.Thermal.HeatTransfer.Interfaces.PartialElementaryConditionalHeatPortWithoutT;
    protected
      Modelica.SIunits.Force f_c "Spring force";
      Modelica.SIunits.Force f_d "Damping force";
    equation
      s_rel = flange_b.s - flange_a.s;
      v_rel = der(s_rel);
      flange_b.f = f;
      flange_a.f = -f;
      f_c = c * (s_rel - s_rel0);
      f_d = d * v_rel;
      f = f_c + f_d;
      lossPower = f_d * v_rel;
      annotation(
        Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Line(origin = {-80, -10}, points = {{0, 70}, {0, -70}}, color = {181, 176, 179}), Polygon(origin = {-79, 70}, lineColor = {211, 209, 209}, fillColor = {216, 216, 216}, fillPattern = FillPattern.Solid, points = {{-1, 10}, {-9, -10}, {9, -10}, {-1, 10}}), Line(origin = {0, 88}, points = {{0, 8}, {0, -8}}), Line(origin = {0, -88}, points = {{0, 8}, {0, -8}}), Line(origin = {0, 80}, points = {{-40, 0}, {40, 0}}), Line(origin = {0, -80}, points = {{-40, 0}, {40, 0}}), Line(origin = {-40, 20}, points = {{-10, 0}, {10, 0}}), Line(origin = {-40, 50}, points = {{0, 30}, {0, -30}}), Line(origin = {-40, 16}, points = {{-12, 14}, {-12, -16}, {12, -16}, {12, 14}}), Line(origin = {-40, -40}, points = {{0, 40}, {0, -40}}), Line(origin = {39.8701, 0}, points = {{0.129947, 80}, {0.129947, 60}, {-19.8701, 40}, {20.1299, 20}, {-19.8701, 0}, {20.1299, -20}, {-19.8701, -40}, {0.129947, -60}, {0.129947, -80}}), Text(origin = {-1, 132}, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-119, -32}, {119, 32}}, textString = "%name")}));
    end DampedSpring_Vertical;

    model SuspendedMass
      //Rigid body with mass and two flanges
      Modelica.SIunits.Position s "Absolute position of center of component (s = flange_a.s + L/2 = flange_b.s - L/2)";
      parameter Modelica.SIunits.Length H(start = 0.2) "Length of component, from left flange to right flange (= flange_b.s - flange_a.s)";
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a "Bottom flange of translational component" annotation(
        Placement(transformation(extent = {{-10, -90}, {10, -70}})));
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Top flange of translational component" annotation(
        Placement(transformation(extent = {{-10, 90}, {10, 70}})));
      constant Real g = Modelica.Constants.g_n;
      parameter Modelica.SIunits.Mass m = 50;
      Modelica.SIunits.Velocity v;
      Modelica.SIunits.Acceleration a;
    equation
      flange_a.s = s - H / 2;
      flange_b.s = s + H / 2;
      v = der(s);
      a = der(v);
      flange_b.f - m * g = m * a;
      annotation(
        Diagram(graphics = {Rectangle(extent = {{-100, 80}, {100, -80}}), Text(origin = {5, 0}, extent = {{-53, 18}, {53, -18}}, textString = "Mass")}),
        Icon(graphics = {Rectangle(extent = {{-100, 80}, {100, -80}}), Text(origin = {-2, 6}, extent = {{-56, 34}, {56, -34}}, textString = "Mass"), Line(origin = {10.0883, -10.0883}, points = {{-70.0883, 50.0883}, {49.9117, 50.0883}, {69.9117, -49.9117}}), Line(origin = {0.00970966, -9.91165}, points = {{-60.0097, 49.9117}, {-80.0097, -50.0883}, {79.9903, -50.0883}}), Line(origin = {0, 50}, points = {{-20, -10}, {-20, 10}, {20, 10}, {20, -10}, {20, -10}}), Text(origin = {-8, 122}, lineColor = {0, 0, 255}, extent = {{-112, -22}, {112, 22}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end SuspendedMass;
  end MiscellaneousComponents;

  package Functions
    function waveNumberIterator
      input Real d "Water depth";
      input Real omega[:] "Wave frequency";
      output Real k[size(omega, 1)] "Wave number";
    protected
      constant Real g = Modelica.Constants.g_n;
      constant Real pi = Modelica.Constants.pi;
      parameter Integer n = size(omega, 1);
      Real T[size(omega, 1)] "Wave period";
      Real L0[size(omega, 1)] "Deepwater wavelength";
      Real Ll(start = 0, fixed = true) "Intermediate loop value";
      Real Llc(start = 0, fixed = true) "Intermediate loop comparator value";
      Real L[size(omega, 1)] "Iterated wave length";
    algorithm
      T := 2 * pi ./ omega;
      L0 := g * T .^ 2 / (2 * pi);
      for i in 1:size(omega, 1) loop
        Ll := L0[i];
        Llc := 0;
        while abs(Llc - Ll) > 0.001 loop
          Llc := Ll;
          L[i] := g * T[i] ^ 2 / (2 * pi) * tanh(2 * pi / Ll * d);
          Ll := L[i];
        end while;
      end for;
      k := 2 * pi ./ L;
    end waveNumberIterator;

    function morisonForceCydlBuoy
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real x_b;
      input Real v_b;
      input Real a_b;
      input Real epsilon[:];
      input Real zeta0i[:];
      input Real rho_w;
      input Real SSE_X;
      input Real r;
      input Real h;
      input Real Cd_b;
      input Real Cma_b;
      input Real z3;
      input Real pz;
      input Real zcg[:];
      input Real Ucg[:];
      output Real Mf;
    protected
      constant Real pi = Modelica.Constants.pi;
      parameter Real npi = ceil(abs(SSE_X - z3) / pz);
      parameter Integer np = integer(npi);
      parameter Real npai = ceil(abs(z3 + h) / pz);
      parameter Integer npa = integer(npai);
      parameter Real Cm_b = 1 + Cma_b;
      parameter Real Ai = rho_w * pi * (2 * r) ^ 2 / 4;
      parameter Real Ad = 0.5 * rho_w * 2 * r;
      Real MCma;
      Real zb[np + 1];
      Real zbm[np];
      Real U_c[np];
      Real u_w[np];
      Real a_w[np];
      Real Mfi[np];
      Real zba[npa + 1];
      Real zbam[npa];
      Real U_ca[npa];
      Real u_wa[npa];
      Real a_wa[npa];
      Real Mfia[npa];
      Real zb1[2];
      Real zb1m[1];
      Real U_c1[1];
      Real u_w1[1];
      Real a_w1[1];
      Real Mf1;
    algorithm
      if sum(zeta0i) > 0 then
        MCma := Cma_b * Ai * a_b;
      else
        MCma := 0;
      end if;
      if z3 >= SSE_X then
        Mf := 0;
      elseif sum(Ucg .^ 2) == 0 and sum(zeta0i) == 0 then
        Mf := -100000 * Cd_b * Ad * abs(v_b) * v_b * abs(z3);
      elseif np > 1 then
        if z3 + h > SSE_X then
          zb[1] := z3;
          for i in 2:np loop
            zb[i] := z3 + (i - 1) * pz;
          end for;
          zb[np + 1] := SSE_X;
          for i in 1:np loop
            zbm[i] := (zb[i] + zb[i + 1]) / 2;
          end for;
          U_c := OceanEngineering.Functions.linearInterpolatorMV(zcg + SSE_X * ones(size(zcg, 1)), Ucg, zbm);
          (u_w, a_w) := waveKinematics(t, d, omega, T, k, x_b, epsilon, zeta0i, SSE_X, zbm);
          for i in 1:np loop
            Mfi[i] := (Cm_b * Ai * a_w[i] - MCma + Cd_b * Ad * abs(u_w[i] + U_c[i] - v_b) * (u_w[i] + U_c[i] - v_b)) * abs(zb[i + 1] - zb[i]);
          end for;
          Mf := sum(Mfi);
        else
          zba[1] := z3;
          for i in 2:npa loop
            zba[i] := z3 + (i - 1) * pz;
            if zba[i] >= z3 + h - pz then
              break;
            end if;
          end for;
          zba[npa + 1] := z3 + h;
          for i in 1:npa loop
            zbam[i] := (zba[i] + zba[i + 1]) / 2;
          end for;
          U_ca := OceanEngineering.Functions.linearInterpolatorMV(zcg + SSE_X * ones(size(zcg, 1)), Ucg, zbam);
          (u_wa, a_wa) := OceanEngineering.Functions.waveKinematics(t, d, omega, T, k, x_b, epsilon, zeta0i, SSE_X, zbam);
          for i in 1:npa loop
            Mfia[i] := (Cm_b * Ai * a_wa[i] - MCma + Cd_b * Ad * abs(u_wa[i] + U_ca[i] - v_b) * (u_wa[i] + U_ca[i] - v_b)) * abs(zba[i + 1] - zba[i]);
          end for;
          Mf := sum(Mfia);
        end if;
      else
        zb1[1] := z3;
        zb1[2] := SSE_X;
        zb1m[1] := (zb1[1] + zb1[2]) / 2;
        U_c1 := OceanEngineering.Functions.linearInterpolatorMV(zcg + SSE_X * ones(size(zcg, 1)), Ucg, zb1m);
        (u_w1, a_w1) := OceanEngineering.Functions.waveKinematics(t, d, omega, T, k, x_b, epsilon, zeta0i, SSE_X, zb1m);
        Mf1 := (Cm_b * Ai * a_w1[1] - MCma + Cd_b * Ad * abs(u_w1[1] + U_c1[1] - v_b) * (u_w1[1] + U_c1[1] - v_b)) * (zb1[2] - zb1[1]);
        Mf := Mf1;
      end if;
    end morisonForceCydlBuoy;

    function linearInterpolatorMV
      input Real x[:];
      input Real y[:];
      input Real p[:];
      output Real q[size(p, 1)];
    protected
      Integer i, j;
    algorithm
      for i in 1:size(p, 1) loop
        if p[i] == x[size(x, 1)] then
          q[i] := y[size(x, 1)];
        else
          j := 1;
          while p[i] >= x[j + 1] loop
            j := j + 1;
          end while;
          q[i] := (y[j + 1] - y[j]) / (x[j + 1] - x[j]) * (p[i] - x[j]) + y[j];
        end if;
      end for;
    end linearInterpolatorMV;

    function waveKinematics
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real x;
      input Real epsilon[:];
      input Real zeta0i[:];
      input Real SSE_X;
      input Real zbm[:] "z coordinates of panel midpoint";
      output Real u[size(zbm, 1)] "effective velocity at zbm";
      output Real a[size(zbm, 1)] "effective acceleration at zbm";
    protected
      constant Real pi = Modelica.Constants.pi;
      Real theta[:];
      Real alpha[:];
      Real H[size(omega, 1)];
      Real ui[size(omega, 1)];
      Real ai[size(omega, 1)];
    algorithm
      H := 2 * zeta0i;
      theta := k * x - omega * t - 2 * pi * epsilon;
      for i in 1:size(zbm, 1) loop
        alpha := cosh(k * (zbm[i] + SSE_X + d)) ./ sinh(k * d);
        ui := pi * H ./ T .* alpha .* cos(theta);
        ai := 2 * pi ^ 2 * H ./ T .^ 2 .* alpha .* sin(theta);
        u[i] := sum(ui);
        a[i] := sum(ai);
      end for;
    end waveKinematics;

    function linearInterpolatorSV
      input Real x[:];
      input Real y[:];
      input Real p;
      output Real q;
    protected
      Integer j;
    algorithm
      if p == x[size(x, 1)] then
        q := y[size(x, 1)];
      else
        j := 1;
        while p >= x[j + 1] loop
          j := j + 1;
        end while;
        q := (y[j + 1] - y[j]) / (x[j + 1] - x[j]) * (p - x[j]) + y[j];
      end if;
    end linearInterpolatorSV;

    function sSE_X_cat
      input Real t;
      input Real omega[:];
      input Real k[:];
      input Real zeta0i[:];
      input Real epsilon[:];
      input Real x_cat_mor;
      output Real SSE_X_cat_mor;
    protected
      constant Real pi = Modelica.Constants.pi;
    algorithm
      SSE_X_cat_mor := sum(zeta0i .* cos(k * x_cat_mor - omega * t - 2 * pi * epsilon));
    end sSE_X_cat;

    function catThIterator
      input Real l;
      input Real x[:];
      input Real h;
      input Real spm_chain_sub;
      output Real Th[size(x, 1)];
    protected
      constant Real g = Modelica.Constants.g_n;
      Real Thl;
      Real Thlc;
      Real ls[size(x, 1)];
      Real X[size(x, 1)];
    algorithm
      for j in 1:size(x, 1) loop
        if x[j] == 0 then
          Th[j] := 0;
        else
          Thl := 500;
          Thlc := 0;
          while abs(Thlc - Thl) > 0.00001 loop
            Thlc := Thl;
            Th[j] := x[j] * spm_chain_sub * g / Modelica.Math.acosh(1 + spm_chain_sub * g * h / Thlc);
            Thl := Th[j];
          end while;
        end if;
      end for;
    end catThIterator;

    function catXIterator
      input Real l;
      input Real h;
      input Real spm_chain_sub;
      input Real Th[:];
      input Real x[:];
      output Real X[size(Th, 1)];
    protected
      constant Real g = Modelica.Constants.g_n;
      Real ls[size(Th, 1)];
    algorithm
      for i in 1:size(Th, 1) loop
        if Th[i] == 0 then
          ls[i] := h;
          X[i] := l - ls[i];
        else
          ls[i] := h * (1 + 2 * (Th[i] / (spm_chain_sub * g * h))) ^ 0.5;
          X[i] := l - ls[i] + x[i];
        end if;
      end for;
    end catXIterator;

    function wave_uCalculator
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real zeta0i[:];
      input Real epsilon[:];
      input Real x;
      input Real z;
      input Real SSE_X;
      output Real u;
    protected
      constant Real pi = Modelica.Constants.pi;
      Real theta[:];
      Real alpha[:];
      Real H[size(omega, 1)];
    algorithm
      H := 2 * zeta0i;
      theta := k * x - omega * t - 2 * pi * epsilon;
      alpha := cosh(k * (z + SSE_X + d)) ./ sinh(k * d);
      u := sum(pi * H ./ T .* alpha .* cos(theta));
    end wave_uCalculator;

    function wave_wCalculator
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real zeta0i[:];
      input Real epsilon[:];
      input Real x;
      input Real z;
      input Real SSE_X;
      output Real w;
    protected
      constant Real pi = Modelica.Constants.pi;
      Real theta[:];
      Real alpha[:];
      Real H[size(omega, 1)];
    algorithm
      H := 2 * zeta0i;
      theta := k * x - omega * t - 2 * pi * epsilon;
      alpha := sinh(k * (z + SSE_X + d)) ./ sinh(k * d);
      w := sum(pi * H ./ T .* alpha .* sin(theta));
    end wave_wCalculator;

    function wave_auCalculator
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real zeta0i[:];
      input Real epsilon[:];
      input Real x;
      input Real z;
      input Real SSE_X;
      output Real au;
    protected
      constant Real pi = Modelica.Constants.pi;
      Real theta[:];
      Real alpha[:];
      Real H[size(omega, 1)];
    algorithm
      H := 2 * zeta0i;
      theta := k * x - omega * t - 2 * pi * epsilon;
      alpha := cosh(k * (z + SSE_X + d)) ./ sinh(k * d);
      au := sum(2 * pi ^ 2 * H ./ T .^ 2 .* alpha .* sin(theta));
    end wave_auCalculator;

    function wave_awCalculator
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real zeta0i[:];
      input Real epsilon[:];
      input Real x;
      input Real z;
      input Real SSE_X;
      output Real aw;
    protected
      constant Real pi = Modelica.Constants.pi;
      Real theta[:];
      Real alpha[:];
      Real H[size(omega, 1)];
    algorithm
      H := 2 * zeta0i;
      theta := k * x - omega * t - 2 * pi * epsilon;
      alpha := sinh(k * (z + SSE_X + d)) ./ sinh(k * d);
      aw := sum((-2 * pi ^ 2 * H ./ T .^ 2) .* alpha .* cos(theta));
    end wave_awCalculator;

    function randomNumberGenerator
      input Integer ls = 614657;
      input Integer gs = 30020;
      input Integer n = 100;
      output Real r64[n];
    protected
      Integer state64[2](each start = 0, each fixed = true);
    algorithm
      state64[1] := 0;
      state64[2] := 0;
      for i in 1:n loop
        if i == 1 then
          state64 := Modelica.Math.Random.Generators.Xorshift64star.initialState(ls, gs);
          r64[i] := 0;
        else
          (r64[i], state64) := Modelica.Math.Random.Generators.Xorshift64star.random(pre(state64));
        end if;
      end for;
    end randomNumberGenerator;

    function frequencySelector
      input Real omega_min;
      input Real omega_max;
      input Real epsilon[:];
      output Real omega[size(epsilon, 1)];
    protected
      parameter Real ref_omega[size(epsilon, 1)] = omega_min:(omega_max - omega_min) / (size(epsilon, 1) - 1):omega_max;
    algorithm
      omega[1] := omega_min;
      for i in 2:size(epsilon, 1) - 1 loop
        omega[i] := ref_omega[i] + epsilon[i] * omega_min;
      end for;
      omega[size(epsilon, 1)] := omega_max;
    end frequencySelector;

    function spectrumGenerator_PM
      input Real Hs = 2 "Significant wave height";
      input Real omega[:] "Frequency components";
      output Real spec[size(omega, 1)] "Spectral values for input frequencies";
    protected
      constant Real pi = Modelica.Constants.pi;
      constant Real g = Modelica.Constants.g_n;
    algorithm
      for i in 1:size(omega, 1) loop
        spec[i] := 0.0081 * g ^ 2 / omega[i] ^ 5 * exp(-0.032 * (g / (Hs * omega[i] ^ 2)) ^ 2);
      end for;
    end spectrumGenerator_PM;

    function morisonForceBoxPlatform
      input Real t;
      input Real d;
      input Real omega[:];
      input Real T[:];
      input Real k[:];
      input Real x;
      input Real v;
      input Real a;
      input Real epsilon[:];
      input Real zeta0i[:];
      input Real rho_w;
      input Real SSE_X;
      input Real L;
      input Real B;
      input Real H;
      input Real Awp;
      input Real Cd;
      input Real Cma;
      input Real z2;
      input Real z3;
      input Real pz;
      input Real zcg[:];
      input Real Ucg[:];
      output Real Mf;
    protected
      parameter Real npi = ceil(abs(SSE_X - z3) / pz);
      parameter Integer np = integer(npi);
      parameter Real npai = ceil(abs(z3 + H) / pz);
      parameter Integer npa = integer(npai);
      parameter Real Cm = 1 + Cma;
      parameter Real Ai = rho_w * Awp;
      parameter Real Ad = 0.5 * rho_w * B;
      Real MCma;
      Real zb[np + 1];
      Real zbm[np];
      Real U_c[np];
      Real u_w[np];
      Real a_w[np];
      Real Mfi[np];
      Real zba[npa + 1];
      Real zbam[npa];
      Real U_ca[npa];
      Real u_wa[npa];
      Real a_wa[npa];
      Real Mfia[npa];
      Real zb1[2];
      Real zb1m[1];
      Real U_c1[1];
      Real u_w1[1];
      Real a_w1[1];
      Real Mf1;
    algorithm
      if sum(zeta0i) > 0 then
        MCma := Cma * Ai * a;
      else
        MCma := 0;
      end if;
      if z3 >= SSE_X then
        Mf := 0;
      elseif sum(Ucg .^ 2) == 0 and sum(zeta0i) == 0 then
        Mf := 0;
      elseif np > 1 then
        if z3 + H > SSE_X then
          zb[1] := z3;
          for i in 2:np loop
            zb[i] := z3 + (i - 1) * pz;
          end for;
          zb[np + 1] := SSE_X;
          for i in 1:np loop
            zbm[i] := (zb[i] + zb[i + 1]) / 2;
          end for;
          U_c := linearInterpolatorMV(zcg + SSE_X * ones(size(zcg, 1)), Ucg, zbm);
          (u_w, a_w) := waveKinematics(t, d, omega, T, k, x, epsilon, zeta0i, SSE_X, zbm);
          for i in 1:np loop
            Mfi[i] := (Cm * Ai * a_w[i] - MCma + Cd * Ad * abs(u_w[i] + U_c[i] - v) * (u_w[i] + U_c[i] - v)) * abs(zb[i + 1] - zb[i]);
          end for;
          Mf := sum(Mfi);
        else
          zba[1] := z3;
          for i in 2:npa loop
            zba[i] := z3 + (i - 1) * pz;
            if zba[i] >= z3 + H - pz then
              break;
            end if;
          end for;
          zba[npa + 1] := z3 + H;
          for i in 1:npa loop
            zbam[i] := (zba[i] + zba[i + 1]) / 2;
          end for;
          U_ca := linearInterpolatorMV(zcg + SSE_X * ones(size(zcg, 1)), Ucg, zbam);
          (u_wa, a_wa) := waveKinematics(t, d, omega, T, k, x, epsilon, zeta0i, SSE_X, zbam);
          for i in 1:npa loop
            Mfia[i] := (Cm * Ai * a_wa[i] - MCma + Cd * Ad * abs(u_wa[i] + U_ca[i] - v) * (u_wa[i] + U_ca[i] - v)) * abs(zba[i + 1] - zba[i]);
          end for;
          Mf := sum(Mfia);
        end if;
      else
        zb1[1] := z3;
        zb1[2] := SSE_X;
        zb1m[1] := (zb1[1] + zb1[2]) / 2;
        U_c1 := linearInterpolatorMV(zcg + SSE_X * ones(size(zcg, 1)), Ucg, zb1m);
        (u_w1, a_w1) := waveKinematics(t, d, omega, T, k, x, epsilon, zeta0i, SSE_X, zb1m);
        Mf1 := (Cm * Ai * a_w1[1] - MCma + Cd * Ad * abs(u_w1[1] + U_c1[1] - v) * (u_w1[1] + U_c1[1] - v)) * (zb1[2] - zb1[1]);
        Mf := Mf1;
      end if;
    end morisonForceBoxPlatform;
  end Functions;

  package Connectors
    connector WaveDataConnector
      parameter Integer n_omega_i = 100 "number of wave components";
      Modelica.Blocks.Interfaces.RealOutput d;
      Modelica.Blocks.Interfaces.RealOutput rho_w;
      Modelica.Blocks.Interfaces.RealOutput omega[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput T[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput k[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput epsilon[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput zeta0i[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput SSE_X0;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end WaveDataConnector;

    expandable connector EnvironmentBus
      annotation(
        Icon(graphics = {Rectangle(fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Line(origin = {0, 80}, points = {{-80, 0}, {80, 0}}), Line(origin = {-0.566026, -0.471698}, points = {{-80, 0}, {80, 0}}), Line(origin = {-0.471698, 39.9057}, points = {{-80, 0}, {80, 0}}), Line(origin = {0.51888, -40}, points = {{-80, 0}, {80, 0}}), Line(origin = {-0.707535, -79.8585}, points = {{-80, 0}, {80, 0}}), Rectangle(origin = {-2, -1}, fillColor = {255, 0, 127}, fillPattern = FillPattern.Solid, extent = {{-52, 91}, {52, -91}}), Rectangle(origin = {-100, 0}, fillColor = {85, 85, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Line(origin = {-86, 50}, points = {{6, 30}, {-6, -30}}), Line(origin = {-83.1673, 29.8327}, points = {{2.66543, 10}, {-3.33457, -10}}), Line(origin = {-88.171, -50}, points = {{7, -30}, {-5, 30}}), Line(origin = {-84, -30}, points = {{4, -10}, {-2, 10}}), Rectangle(origin = {100, 0}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}}), Line(origin = {87.8284, 50}, points = {{-7, 30}, {7, -30}}), Line(origin = {83, 30}, points = {{-3, 10}, {3, -10}}), Line(origin = {84, -30}, points = {{-4, -10}, {4, 10}}), Line(origin = {87.8284, -50}, points = {{-9, -30}, {7, 30}}), Text(origin = {1, 113}, lineColor = {0, 0, 255}, extent = {{-150, 60}, {150, -10}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end EnvironmentBus;

    connector EnvironmentBuoyDataConnector
      parameter Integer n_omega_i = 100 "number of wave components";
      Modelica.Blocks.Interfaces.RealInput xinit;
      Modelica.Blocks.Interfaces.RealInput omega[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput T[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput k[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput epsilon[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput zeta0i[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput zcg[4];
      Modelica.Blocks.Interfaces.RealInput Ucg[4];
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end EnvironmentBuoyDataConnector;

    connector EnvironmentMooringDataConnector
      parameter Integer n_omega_i = 100 "number of wave components";
      Modelica.Blocks.Interfaces.RealOutput xinit;
      Modelica.Blocks.Interfaces.RealInput omega[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput T[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput k[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput epsilon[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput zeta0i[n_omega_i];
      Modelica.Blocks.Interfaces.RealInput zcg[4];
      Modelica.Blocks.Interfaces.RealInput Ucg[4];
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end EnvironmentMooringDataConnector;

    connector CurrentDataConnector
      Modelica.Blocks.Interfaces.RealOutput zcg[4];
      Modelica.Blocks.Interfaces.RealOutput Ucg[4];
      annotation(
        Icon(graphics = {Rectangle(fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end CurrentDataConnector;
  end Connectors;

  package SampleSimulations
    extends Modelica.Icons.ExamplesPackage;

    package FreeFloatingCylinder
      connector EnvironmentBuoyDataConnector_Free
        parameter Integer n_omega_i = 100 "number of wave components";
        Modelica.Blocks.Interfaces.RealInput omega[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput T[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput k[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput epsilon[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput zeta0i[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput zcg[4];
        Modelica.Blocks.Interfaces.RealInput Ucg[4];
        annotation(
          Icon(graphics = {Rectangle(fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
      end EnvironmentBuoyDataConnector_Free;

      model CylindricalBuoy_Free
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_a TopHook(s = z2) "Hook for spring" annotation(
          Placement(transformation(extent = {{-20, 120}, {20, 80}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a Fairlead[2](s = {x, z3}) annotation(
          Placement(transformation(extent = {{-20, -120}, {20, -80}})));
        OceanEngineering.SampleSimulations.FreeFloatingCylinder.EnvironmentBuoyDataConnector_Free ebdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{-80, -20}, {-120, 20}}, rotation = 0), iconTransformation(extent = {{-118, -20}, {-78, 20}}, rotation = 0)));
        constant Real g = Modelica.Constants.g_n "Value of acceleration due to gravity";
        constant Real pi = Modelica.Constants.pi "Value of pi";
        parameter Integer nOmega = 1 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        parameter Real d = 100;
        parameter Real rho_w = 1025;
        parameter Modelica.SIunits.Length r = 0.6 "Radius of the buoy";
        parameter Modelica.SIunits.Length h = 2 "Height of the buoy";
        parameter Modelica.SIunits.Mass m_s = 350 "Structural mass of the buoy";
        parameter Modelica.SIunits.Mass m_b = 500 "Ballast mass";
        parameter Modelica.SIunits.Length pz = 2 "Depth of panel to calculate current force & Morison Force";
        parameter Modelica.SIunits.Length zKG = 1.25 "KG of buoy";
        parameter Real Cma_x = 1 "Added mass coefficient";
        parameter Real Cma_z = 1 "Added mass coefficient";
        parameter Real Cd_x = 1 "Drag coefficient of buoy";
        parameter Real Cd_z = 1 "Drag coefficient of buoy";
        parameter Real K_x(unit = "N m^-1") = 0 "Coefficient of Stiffness";
        parameter Real K_z(unit = "N m^-1") = rho_w * g * pi * r ^ 2 "Coefficient of Stiffness";
        parameter Real C_x(unit = "kg/s") = 0 "Damping in x";
        parameter Real C_z(unit = "kg/s") = 3100 "Damping in z";
        parameter Modelica.SIunits.Mass M = m_s + m_b "Dry mass of the buoy";
        parameter Modelica.SIunits.Area Awp = pi * r ^ 2 "Waterplane area of buoy";
        parameter Integer mor_flg = 0 " '1' if moored, '0' if free floating";
        parameter Modelica.SIunits.Length z_s = M / (Awp * rho_w) "Static draught";
        parameter Modelica.SIunits.Length z_fb = zKG - z_s "VCG of body wrt SWL";
        parameter Modelica.SIunits.Mass spm_chain = 14 "Specific mass of the mooring line";
        parameter Modelica.SIunits.Density rho_mat_chain = 7800 "Density of mooring line material";
        parameter Modelica.SIunits.Mass spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Modelica.SIunits.Length z_m = z_s+mor_flg*(spm_chain_sub * (d-z_s) * g / K_z) "Approximate change in static draught when mooring chain is attached";
        parameter Modelica.SIunits.Length z_cb = -z_m / 2 "Static centre of buoyancy";
        Modelica.SIunits.Length SSE_X "Sea surface elevation corresponding to buoy CG X coordinate";
        Modelica.SIunits.Acceleration a_w "Water particle acceleration in z direction";
        Modelica.SIunits.Length x "x co-ordinate of VCG";
        Modelica.SIunits.Length z "Displacement of body VCG from z_fb";
        Modelica.SIunits.Velocity v_x "Velocity of body VCG in horizontal direction";
        Modelica.SIunits.Acceleration a_x "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Velocity v_z "Velocity of body VCG in vertical direction";
        Modelica.SIunits.Acceleration a_z "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Length z1 "Position of body VCG";
        Modelica.SIunits.Length z2 "Position of deck";
        Modelica.SIunits.Length z3 "Position of keel";
        Modelica.SIunits.Force MF_x "Morison force due to current and waves";
      initial equation
        z3=-z_m;
        der(z) = 0;
        x = 0;
        der(x) = 0;
      equation
        SSE_X = sum(ebdc.zeta0i .* cos(ebdc.k * x - ebdc.omega * time - 2 * pi * ebdc.epsilon));
        a_w = OceanEngineering.Functions.wave_awCalculator(time, d, ebdc.omega, ebdc.T, ebdc.k, ebdc.zeta0i, ebdc.epsilon, x, z_cb, SSE_X);
        v_z = der(z);
        a_z = der(v_z);
        M * a_z + C_z * (SSE_X - z3) * v_z + K_z * z = rho_w * g * Awp * SSE_X + Cma_x * M * a_w + TopHook.f + Fairlead[2].f;
        z1 = z_fb + z;
        z2 = z3 + h;
        z3 = -z_s + z;
        MF_x = OceanEngineering.Functions.morisonForceCydlBuoy(time, d, ebdc.omega, ebdc.T, ebdc.k, x, v_x, a_x, ebdc.epsilon, ebdc.zeta0i, rho_w, SSE_X, r, h, Cd_x, Cma_x, z3, pz, ebdc.zcg, ebdc.Ucg);
        v_x = der(x);
        a_x = der(v_x);
        M * a_x + C_x * M * (SSE_X - z3) * v_x + K_x * x = Fairlead[1].f + MF_x;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -50}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 50}, {100, -50}}), Rectangle(origin = {0, 3}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 53}, {28, -53}}), Ellipse(origin = {0, 57}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360), Ellipse(origin = {0, -51}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 300, Tolerance = 1e-05, Interval = 0.5));
      end CylindricalBuoy_Free;






      model FreeFloatingCylinderInRegularWaves
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(Hr = 1, Tr = 7) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Free cylindricalBuoy_Free1 annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(Uf = {0, 0, 0, 0}) annotation(
          Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-40, 30}, {-27, 30}, {-27, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-40, -30}, {-20, -30}, {-20, 0}, {0, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T) annotation(
          Line);
        connect(cylindricalBuoy_Free1.ebdc.omega, environmentBus1.omega) annotation(
          Line(points = {{0, 0}, {30, 0}, {30, 0}, {30, 0}}, thickness = 0.5));
        connect(cylindricalBuoy_Free1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy_Free1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy_Free1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy_Free1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy_Free1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy_Free1.ebdc.Ucg, environmentBus1.Ucg);
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.1));
      end FreeFloatingCylinderInRegularWaves;

      model FreeFloatingCylinderInRegularWavesAndCurrent
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(Hr = 1, Tr = 7,d=50,Tdel=0) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Free cylindricalBuoy_Free1(d=50) annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp=20) annotation(
          Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-40, -30}, {-20, -30}, {-20, 0}, {0, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-40, 30}, {-15, 30}, {-15, 0}, {0, 0}}));
        connect(cylindricalBuoy_Free1.ebdc.omega, environmentBus1.omega) annotation(
          Line(points = {{0, 0}, {30, 0}, {30, 0}, {30, 0}}, thickness = 0.5));
        connect(cylindricalBuoy_Free1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy_Free1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy_Free1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy_Free1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy_Free1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy_Free1.ebdc.Ucg, environmentBus1.Ucg);
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.1));
      end FreeFloatingCylinderInRegularWavesAndCurrent;



      model FreeFloatingCylinderInCurrent
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(Hr = 0, Tr = 7, d = 50, Tdel = 0) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Free cylindricalBuoy_Free1(d = 50) annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-40, -30}, {-20, -30}, {-20, 0}, {0, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-40, 30}, {-15, 30}, {-15, 0}, {0, 0}}));
        connect(cylindricalBuoy_Free1.ebdc.omega, environmentBus1.omega) annotation(
          Line(points = {{0, 0}, {30, 0}, {30, 0}, {30, 0}}, thickness = 0.5));
        connect(cylindricalBuoy_Free1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy_Free1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy_Free1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy_Free1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy_Free1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy_Free1.ebdc.Ucg, environmentBus1.Ucg);
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.1));
      end FreeFloatingCylinderInCurrent;




    end FreeFloatingCylinder;

    package FixedCylinder
      model CylindricalBuoy_Fixed
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_a TopHook(s = z2) "Hook for spring" annotation(
          Placement(transformation(extent = {{-20, 120}, {20, 80}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a Fairlead[2](s = {x, z3}) annotation(
          Placement(transformation(extent = {{-20, -120}, {20, -80}})));
        OceanEngineering.SampleSimulations.FixedCylinder.EnvironmentBuoyDataConnector_Fixed ebdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{-80, -20}, {-120, 20}}, rotation = 0), iconTransformation(extent = {{-118, -20}, {-78, 20}}, rotation = 0)));
        constant Real g = Modelica.Constants.g_n "Value of acceleration due to gravity";
        constant Real pi = Modelica.Constants.pi "Value of pi";
        parameter Integer nOmega = 1 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        parameter Real d = 100;
        parameter Real rho_w = 1025;
        parameter Modelica.SIunits.Length r = 0.6 "Radius of the buoy";
        parameter Modelica.SIunits.Length h = 2 "Height of the buoy";
        parameter Modelica.SIunits.Mass m_s = 350 "Structural mass of the buoy";
        parameter Modelica.SIunits.Mass m_b = 500 "Ballast mass";
        parameter Modelica.SIunits.Length pz = 2 "Depth of panel to calculate current force & Morison Force";
        parameter Modelica.SIunits.Length zKG = 1.25 "KG of buoy";
        parameter Real Cma_x = 1 "Added mass coefficient";
        parameter Real Cma_z = 1 "Added mass coefficient";
        parameter Real Cd_x = 1 "Drag coefficient of buoy";
        parameter Real Cd_z = 1 "Drag coefficient of buoy";
        parameter Real K_x(unit = "N m^-1") = 0 "Coefficient of Stiffness";
        parameter Real K_z(unit = "N m^-1") = rho_w * g * pi * r ^ 2 "Coefficient of Stiffness";
        parameter Real C_x(unit = "kg/s") = 0 "Damping in x";
        parameter Real C_z(unit = "kg/s") = 3100 "Damping in z";
        parameter Modelica.SIunits.Mass M = m_s + m_b "Dry mass of the buoy";
        parameter Modelica.SIunits.Area Awp = pi * r ^ 2 "Waterplane area of buoy";
        parameter Modelica.SIunits.Length z_s = M / (Awp * rho_w) "Static draught";
        parameter Modelica.SIunits.Length z_cg = -z_s / 2 "Static centre of buoyancy";
        parameter Modelica.SIunits.Length z_fb = zKG - z_s "VCG of body wrt SWL";
        parameter Integer mor_flg = 1 " '1' if moored, '0' if free floating";
        parameter Modelica.SIunits.Mass spm_chain = 10 "Specific mass of the mooring line";
        parameter Modelica.SIunits.Density rho_mat_chain = 7800 "Density of mooring line material";
        parameter Modelica.SIunits.Mass spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Modelica.SIunits.Length z_m = spm_chain_sub * (d - z_s) * g / K_z "Approximate change in static draught when mooring chain is attached";
        Modelica.SIunits.Length SSE_X "Sea surface elevation corresponding to buoy CG X coordinate";
        Modelica.SIunits.Acceleration a_w "Water particle acceleration in z direction";
        Modelica.SIunits.Length x "x co-ordinate of VCG";
        Modelica.SIunits.Length z "Displacement of body VCG from z_fb";
        Modelica.SIunits.Velocity v_x "Velocity of body VCG in horizontal direction";
        Modelica.SIunits.Acceleration a_x "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Velocity v_z "Velocity of body VCG in vertical direction";
        Modelica.SIunits.Acceleration a_z "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Length z1 "Position of body VCG";
        Modelica.SIunits.Length z2 "Position of deck";
        Modelica.SIunits.Length z3 "Position of keel";
        Modelica.SIunits.Force MF_x "Morison force due to current and waves";
      equation
        SSE_X = sum(ebdc.zeta0i .* cos(ebdc.k * x - ebdc.omega * time - 2 * pi * ebdc.epsilon));
        a_w = 0;
        v_z = der(z);
        a_z = der(v_z);
        z = 0;
        z1 = z_fb + z;
        z2 = z3 + h;
        z3 = -1;
        MF_x = OceanEngineering.Functions.morisonForceCydlBuoy(time, d, ebdc.omega, ebdc.T, ebdc.k, x, v_x, a_x, ebdc.epsilon, ebdc.zeta0i, rho_w,SSE_X, r, h, Cd_x, Cma_x, z3, pz, ebdc.zcg, ebdc.Ucg);
        v_x = der(x);
        a_x = der(v_x);
        x = 0;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -50}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 50}, {100, -50}}), Rectangle(origin = {0, 3}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 53}, {28, -53}}), Ellipse(origin = {0, 57}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360), Ellipse(origin = {0, -51}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-05, Interval = 0.1));
      end CylindricalBuoy_Fixed;









      connector EnvironmentBuoyDataConnector_Fixed
        parameter Integer n_omega_i = 100 "number of wave components";
        Modelica.Blocks.Interfaces.RealInput omega[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput T[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput k[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput epsilon[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput zeta0i[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput zcg[4];
        Modelica.Blocks.Interfaces.RealInput Ucg[4];
        annotation(
          Icon(graphics = {Rectangle(fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
      end EnvironmentBuoyDataConnector_Fixed;

      model FixedCylinderInCurrent
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(Hr = 0, Tr = 7, d = 50, Tdel = 0) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Fixed cylindricalBuoy_Fixed1(d = 50) annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-40, -30}, {-20, -30}, {-20, 0}, {0, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-40, 30}, {-15, 30}, {-15, 0}, {0, 0}}));
        connect(cylindricalBuoy_Fixed1.ebdc.omega, environmentBus1.omega) annotation(
          Line(points = {{0, 0}, {30, 0}, {30, 0}, {30, 0}}, thickness = 0.5));
        connect(cylindricalBuoy_Fixed1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy_Fixed1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy_Fixed1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy_Fixed1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy_Fixed1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy_Fixed1.ebdc.Ucg, environmentBus1.Ucg);
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.1));
      end FixedCylinderInCurrent;

      model FixedCylinderInRegularWaves
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(Hr = 1, Tr = 7, d = 50, Tdel = 0) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Fixed cylindricalBuoy_Fixed1(d = 50) annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 0, 0}, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-40, -30}, {-20, -30}, {-20, 0}, {0, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-40, 30}, {-15, 30}, {-15, 0}, {0, 0}}));
        connect(cylindricalBuoy_Fixed1.ebdc.omega, environmentBus1.omega) annotation(
          Line(points = {{0, 0}, {30, 0}, {30, 0}, {30, 0}}, thickness = 0.5));
        connect(cylindricalBuoy_Fixed1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy_Fixed1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy_Fixed1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy_Fixed1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy_Fixed1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy_Fixed1.ebdc.Ucg, environmentBus1.Ucg);
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.1));
      end FixedCylinderInRegularWaves;

      model FixedCylinderInRegularWavesAndCurrent
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(Hr = 1, Tr = 7, d = 50, Tdel = 0) annotation(
          Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Fixed cylindricalBuoy_Fixed1(d = 50) annotation(
          Placement(visible = true, transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-40, -30}, {-20, -30}, {-20, 0}, {0, 0}, {0, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T) annotation(
          Line);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-40, 30}, {-15, 30}, {-15, 0}, {0, 0}}));
        connect(cylindricalBuoy_Fixed1.ebdc.omega, environmentBus1.omega) annotation(
          Line(points = {{0, 0}, {30, 0}, {30, 0}, {30, 0}}, thickness = 0.5));
        connect(cylindricalBuoy_Fixed1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy_Fixed1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy_Fixed1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy_Fixed1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy_Fixed1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy_Fixed1.ebdc.Ucg, environmentBus1.Ucg);
        annotation(
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.1));
      end FixedCylinderInRegularWavesAndCurrent;







    end FixedCylinder;

    package MooredCylinder
      model MooredCylinderInCurrentMfC
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 0, Tr = 7, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Floaters.CylindricalBuoy cylindricalBuoy1(d = 50, nOmega = 1, m_b = 500, spm_chain = 10) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_MfC catenaryMooring_MfC1(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.022, spm_chain = 10,Cdn=0.75, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_MfC1.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_MfC1.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_MfC1.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_MfC1.emdc.T, environmentBus1.T);
        connect(catenaryMooring_MfC1.emdc.k, environmentBus1.k);
        connect(catenaryMooring_MfC1.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_MfC1.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_MfC1.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_MfC1.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.5));
      end MooredCylinderInCurrentMfC;







      model MooredCylinderInCurrentMf0
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 0, Tr = 7, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Floaters.CylindricalBuoy cylindricalBuoy1(d = 50, nOmega = 1, m_b = 500, spm_chain = 10) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_Mf0 catenaryMooring_Mf01(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.022, spm_chain =10, Trmp = 20) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_Mf01.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_Mf01.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_Mf01.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_Mf01.emdc.T, environmentBus1.T);
        connect(catenaryMooring_Mf01.emdc.k, environmentBus1.k);
        connect(catenaryMooring_Mf01.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_Mf01.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_Mf01.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_Mf01.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.5));
      end MooredCylinderInCurrentMf0;





      model MooredCylinderInCurrentMfCW
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 0, Tr = 7, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Floaters.CylindricalBuoy cylindricalBuoy1(d = 50, nOmega = 1, m_b = 750, spm_chain = 16) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_MfCW catenaryMooring_MfCW1(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.028, spm_chain = 16, spm_chain_sub = 14, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_MfCW1.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_MfCW1.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_MfCW1.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_MfCW1.emdc.T, environmentBus1.T);
        connect(catenaryMooring_MfCW1.emdc.k, environmentBus1.k);
        connect(catenaryMooring_MfCW1.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_MfCW1.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_MfCW1.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_MfCW1.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.5));
      end MooredCylinderInCurrentMfCW;

      block ThCalculator
        parameter Real d= 50 "water depth";
        parameter Real rho_w=1025 "water density";
        parameter Real L_chain=100 "length of mooring chain";
        parameter Real spm_chain = 43.5 "specific mass of chain in kg/m in air";
        parameter Real d_chain = 0.047 "wire diameter";
        parameter Real rho_mat_chain = 7800 "Density of chain material";
        parameter Real spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged mass per meter of chain";
        parameter Real xmax = sqrt(L_chain ^ 2 - d ^ 2)"maximum displacement of top end";
        parameter Real X0 = L_chain - d"minimum displacement of top end";
        parameter Real x[:] = 0:1:80;
        parameter Real Th[:] = OceanEngineering.Functions.catThIterator(L_chain, x, d, spm_chain_sub);
        parameter Real X[:] = OceanEngineering.Functions.catXIterator(L_chain, d, spm_chain_sub, Th, x);
  annotation(
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 1));
      end ThCalculator;

      model CylindricalBuoy_Fixed
        extends Modelica.Blocks.Icons.Block;
        Modelica.Mechanics.Translational.Interfaces.Flange_a TopHook(s = z2) "Hook for spring" annotation(
          Placement(transformation(extent = {{-20, 120}, {20, 80}})));
        Modelica.Mechanics.Translational.Interfaces.Flange_a Fairlead[2](s = {x, z3}) annotation(
          Placement(transformation(extent = {{-20, -120}, {20, -80}})));
        OceanEngineering.Connectors.EnvironmentBuoyDataConnector ebdc(n_omega_i = nOmega) annotation(
          Placement(visible = true, transformation(extent = {{-80, -20}, {-120, 20}}, rotation = 0), iconTransformation(extent = {{-118, -20}, {-78, 20}}, rotation = 0)));
        constant Real g = Modelica.Constants.g_n "Value of acceleration due to gravity";
        constant Real pi = Modelica.Constants.pi "Value of pi";
        parameter Integer nOmega = 1 "Number of wave components; 1 for regular wave and more than 1 for irregular wave";
        parameter Real d = 100;
        parameter Real rho_w = 1025;
        parameter Modelica.SIunits.Length r = 0.6 "Radius of the buoy";
        parameter Modelica.SIunits.Length h = 2 "Height of the buoy";
        parameter Modelica.SIunits.Mass m_s = 350 "Structural mass of the buoy";
        parameter Modelica.SIunits.Mass m_b = 500 "Ballast mass";
        parameter Modelica.SIunits.Length pz = 2 "Depth of panel to calculate current force & Morison Force";
        parameter Modelica.SIunits.Length zKG = 1.25 "KG of buoy";
        parameter Real Cma_x = 1 "Added mass coefficient";
        parameter Real Cma_z = 1 "Added mass coefficient";
        parameter Real Cd_x = 1 "Drag coefficient of buoy";
        parameter Real Cd_z = 1 "Drag coefficient of buoy";
        parameter Real K_x(unit = "N m^-1") = 0 "Coefficient of Stiffness";
        parameter Real K_z(unit = "N m^-1") = rho_w * g * pi * r ^ 2 "Coefficient of Stiffness";
        parameter Real C_x(unit = "kg/s") = 0 "Damping in x";
        parameter Real C_z(unit = "kg/s") = 3100 "Damping in z";
        parameter Modelica.SIunits.Mass M = m_s + m_b "Dry mass of the buoy";
        parameter Modelica.SIunits.Area Awp = pi * r ^ 2 "Waterplane area of buoy";
        parameter Modelica.SIunits.Length z_s = M / (Awp * rho_w) "Static draught";
        parameter Modelica.SIunits.Length z_cg = -z_s / 2 "Static centre of buoyancy";
        parameter Modelica.SIunits.Length z_fb = zKG - z_s "VCG of body wrt SWL";
        parameter Integer mor_flg = 1 " '1' if moored, '0' if free floating";
        parameter Modelica.SIunits.Mass spm_chain = 10 "Specific mass of the mooring line";
        parameter Modelica.SIunits.Density rho_mat_chain = 7800 "Density of mooring line material";
        parameter Modelica.SIunits.Mass spm_chain_sub = spm_chain - spm_chain / rho_mat_chain * rho_w "Submerged weight per meter of chain";
        parameter Modelica.SIunits.Length z_m = spm_chain_sub * (d - z_s) * g / K_z "Approximate change in static draught when mooring chain is attached";
        Modelica.SIunits.Length SSE_X "Sea surface elevation corresponding to buoy CG X coordinate";
        Modelica.SIunits.Acceleration a_w "Water particle acceleration in z direction";
        Modelica.SIunits.Length x "x co-ordinate of VCG";
        Modelica.SIunits.Length z "Displacement of body VCG from z_fb";
        Modelica.SIunits.Velocity v_x "Velocity of body VCG in horizontal direction";
        Modelica.SIunits.Acceleration a_x "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Velocity v_z "Velocity of body VCG in vertical direction";
        Modelica.SIunits.Acceleration a_z "Acceleration of body VCG in horizontal direction";
        Modelica.SIunits.Length z1 "Position of body VCG";
        Modelica.SIunits.Length z2 "Position of deck";
        Modelica.SIunits.Length z3 "Position of keel";
        Modelica.SIunits.Force MF_x "Morison force due to current and waves";
      equation
        SSE_X = sum(ebdc.zeta0i .* cos(ebdc.k * x - ebdc.omega * time - 2 * pi * ebdc.epsilon));
        a_w = 0;
        v_z = der(z);
        a_z = der(v_z);
        z = 0;
        z1 = z_fb + z;
        z2 = z3 + h;
        z3 = 0;
        MF_x = OceanEngineering.Functions.morisonForceCydlBuoy(time, d, ebdc.omega, ebdc.T, ebdc.k, x, v_x, a_x, ebdc.epsilon, ebdc.zeta0i, rho_w, SSE_X, r, h, Cd_x, Cma_x, z3, pz, ebdc.zcg, ebdc.Ucg);
        v_x = der(x);
        a_x = der(v_x);
        x = 80;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -50}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 50}, {100, -50}}), Rectangle(origin = {0, 3}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 53}, {28, -53}}), Ellipse(origin = {0, 57}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360), Ellipse(origin = {0, -51}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-28, 5}, {28, -5}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-05, Interval = 0.1));
      end CylindricalBuoy_Fixed;






      model CatenaryShape
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 0, Tr = 7, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 0, 0}, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        CylindricalBuoy_Fixed cylindricalBuoy1(d = 50, nOmega = 1, m_b = 500, spm_chain = 10) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_MfC catenaryMooring_MfC1(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.039, spm_chain = 10, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_MfC1.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_MfC1.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_MfC1.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_MfC1.emdc.T, environmentBus1.T);
        connect(catenaryMooring_MfC1.emdc.k, environmentBus1.k);
        connect(catenaryMooring_MfC1.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_MfC1.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_MfC1.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_MfC1.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.5));
      end CatenaryShape;



      model MooredCylinderInWavesAndCurrentMf0
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 1, Tr = 7, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Floaters.CylindricalBuoy cylindricalBuoy1(d = 50, nOmega = 1, m_b = 500, spm_chain = 10) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_Mf0 catenaryMooring_Mf01(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.022, spm_chain = 10, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_Mf01.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_Mf01.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_Mf01.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_Mf01.emdc.T, environmentBus1.T);
        connect(catenaryMooring_Mf01.emdc.k, environmentBus1.k);
        connect(catenaryMooring_Mf01.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_Mf01.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_Mf01.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_Mf01.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.5));
      end MooredCylinderInWavesAndCurrentMf0;



      model MooredCylinderInWavesAndCurrentMfC
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 1, Tr = 7, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 1, 1}, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Floaters.CylindricalBuoy cylindricalBuoy1(d = 50, nOmega = 1, m_b = 500, spm_chain = 10) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_MfC catenaryMooring_MfC1(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.022, spm_chain = 10,Cdn=0.75, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_MfC1.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_MfC1.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_MfC1.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_MfC1.emdc.T, environmentBus1.T);
        connect(catenaryMooring_MfC1.emdc.k, environmentBus1.k);
        connect(catenaryMooring_MfC1.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_MfC1.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_MfC1.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_MfC1.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.5));
      end MooredCylinderInWavesAndCurrentMfC;



      model MooredCylinderInWavesAndCurrentMfCW
        OceanEngineering.Components.Waves.RegularWave.Regular_Airy_Wave regular_Airy_Wave1(d = 50, Hr = 1, Tr = 7, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.Components.CurrentProfiles.CurrentProfile_4pt currentProfile_4pt1(zcg = {-50, -25, -10, 0}, Uf = {0, 0, 0.5, 0.5}, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Floaters.CylindricalBuoy cylindricalBuoy1(d = 50, nOmega = 1, m_b = 500, spm_chain = 10) annotation(
          Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Components.Moorings.CatenaryMooring_MfCW catenaryMooring_MfCW1(nOmega = 1, d = 50, n_lnk = 50, lnk_l = 2, d_chain = 0.022, spm_chain = 10, Trmp = 50) annotation(
          Placement(visible = true, transformation(origin = {30, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(catenaryMooring_MfCW1.shackle, cylindricalBuoy1.Fairlead) annotation(
          Line(points = {{30, -20}, {30, -20}, {30, 20}, {30, 20}}, color = {0, 127, 0}, thickness = 0.5));
        connect(catenaryMooring_MfCW1.emdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, -30}, {0, -30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(catenaryMooring_MfCW1.emdc.omega, environmentBus1.omega);
        connect(catenaryMooring_MfCW1.emdc.T, environmentBus1.T);
        connect(catenaryMooring_MfCW1.emdc.k, environmentBus1.k);
        connect(catenaryMooring_MfCW1.emdc.epsilon, environmentBus1.epsilon);
        connect(catenaryMooring_MfCW1.emdc.zeta0i, environmentBus1.zeta0i);
        connect(catenaryMooring_MfCW1.emdc.zcg, environmentBus1.zcg);
        connect(catenaryMooring_MfCW1.emdc.Ucg, environmentBus1.Ucg);
        connect(cylindricalBuoy1.ebdc.xinit, environmentBus1.xinit) annotation(
          Line(points = {{20, 30}, {0, 30}, {0, 0}, {-20, 0}, {-20, 0}}));
        connect(cylindricalBuoy1.ebdc.omega, environmentBus1.omega);
        connect(cylindricalBuoy1.ebdc.T, environmentBus1.T);
        connect(cylindricalBuoy1.ebdc.k, environmentBus1.k);
        connect(cylindricalBuoy1.ebdc.epsilon, environmentBus1.epsilon);
        connect(cylindricalBuoy1.ebdc.zeta0i, environmentBus1.zeta0i);
        connect(cylindricalBuoy1.ebdc.zcg, environmentBus1.zcg);
        connect(cylindricalBuoy1.ebdc.Ucg, environmentBus1.Ucg);
        connect(currentProfile_4pt1.cdc.zcg, environmentBus1.zcg) annotation(
          Line(points = {{-60, -30}, {-40, -30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(currentProfile_4pt1.cdc.Ucg, environmentBus1.Ucg);
        connect(regular_Airy_Wave1.wdoc.omega, environmentBus1.omega) annotation(
          Line(points = {{-60, 30}, {-40, 30}, {-40, 0}, {-20, 0}, {-20, 0}}));
        connect(regular_Airy_Wave1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.25));
      end MooredCylinderInWavesAndCurrentMfCW;































    end MooredCylinder;
  end SampleSimulations;
  annotation(
    Icon(graphics = {Polygon(origin = {0, 31}, lineColor = {0, 170, 255}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0, points = {{-140, -59}, {-36, -37}, {4, 9}, {54, 13}, {82, -19}, {36, 3}, {16, -53}, {150, -51}, {-140, -59}}, smooth = Smooth.Bezier), Rectangle(origin = {0, -90}, lineColor = {0, 170, 255}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-100, 70}, {100, -10}})}));
end OceanEngineering;
