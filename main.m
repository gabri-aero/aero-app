classdef app1 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        Menu                          matlab.ui.container.Menu
        BasicconfigurationsMenu       matlab.ui.container.Menu
        HalfRankineOvalMenu           matlab.ui.container.Menu
        RankineOvalMenu               matlab.ui.container.Menu
        NonliftingflowoveracylinderMenu  matlab.ui.container.Menu
        LiftingflowoveracylinderMenu  matlab.ui.container.Menu
        BathtubevortexMenu            matlab.ui.container.Menu
        GroundeffectMenu              matlab.ui.container.Menu
        FiniterowofvorticesMenu       matlab.ui.container.Menu
        TabGroup                      matlab.ui.container.TabGroup
        IntroTab                      matlab.ui.container.Tab
        Date05102021Label             matlab.ui.control.Label
        Authors                       matlab.ui.control.Label
        AERODYNAMICFLOWCONFIGURATIONLabel  matlab.ui.control.Label
        Image6                        matlab.ui.control.Image
        Image5                        matlab.ui.control.Image
        BasicflowanalysisTab          matlab.ui.container.Tab
        GridLayout                    matlab.ui.container.GridLayout
        PressurefieldCheckBox         matlab.ui.control.CheckBox
        Panel                         matlab.ui.container.Panel
        GridLayout2                   matlab.ui.container.GridLayout
        IntensitySliderLabel          matlab.ui.control.Label
        IntensitySlider               matlab.ui.control.Slider
        VortexnumberSpinner           matlab.ui.control.Spinner
        VortexnumberSpinnerLabel      matlab.ui.control.Label
        CloseButton                   matlab.ui.control.Button
        SeparationEditField           matlab.ui.control.NumericEditField
        SeparationEditFieldLabel      matlab.ui.control.Label
        AcceptButton                  matlab.ui.control.Button
        UITable                       matlab.ui.control.Table
        AddnewflowButton              matlab.ui.control.Button
        ResetflowsButton              matlab.ui.control.Button
        StreamplotCheckBox            matlab.ui.control.CheckBox
        StreamfunctionCheckBox        matlab.ui.control.CheckBox
        UIAxes                        matlab.ui.control.UIAxes
        JoukowskitransformationTab    matlab.ui.container.Tab
        GridLayout3                   matlab.ui.container.GridLayout
        GridLayout6                   matlab.ui.container.GridLayout
        UIAxes3                       matlab.ui.control.UIAxes
        UIAxes2                       matlab.ui.control.UIAxes
        GridLayout4                   matlab.ui.container.GridLayout
        GridLayout5                   matlab.ui.container.GridLayout
        GridLayout7                   matlab.ui.container.GridLayout
        DLabel                        matlab.ui.control.Label
        LLabel                        matlab.ui.control.Label
        NmLabel_2                     matlab.ui.control.Label
        DragLabel                     matlab.ui.control.Label
        NmLabel                       matlab.ui.control.Label
        LiftLabel                     matlab.ui.control.Label
        TollSliderLabel               matlab.ui.control.Label
        TollSlider                    matlab.ui.control.Slider
        VelocitySlider                matlab.ui.control.Slider
        VelocitySliderLabel           matlab.ui.control.Label
        AoASlider                     matlab.ui.control.Slider
        AoASliderLabel                matlab.ui.control.Label
        GridLayout4_2                 matlab.ui.container.GridLayout
        PressurefieldCheckBox_2       matlab.ui.control.CheckBox
        Image7                        matlab.ui.control.Image
        REditField                    matlab.ui.control.NumericEditField
        REditFieldLabel               matlab.ui.control.Label
        YEditField                    matlab.ui.control.NumericEditField
        YEditFieldLabel               matlab.ui.control.Label
        XEditField                    matlab.ui.control.NumericEditField
        XEditFieldLabel               matlab.ui.control.Label
        PanelmethodanalysisTab        matlab.ui.container.Tab
        GridLayout10                  matlab.ui.container.GridLayout
        NACAairfoilEditField          matlab.ui.control.EditField
        NACAairfoilEditFieldLabel     matlab.ui.control.Label
        PressurefieldCheckBox_3       matlab.ui.control.CheckBox
        ComputepropertiesButton       matlab.ui.control.Button
        StreamplotCheckBox_2          matlab.ui.control.CheckBox
        StreamfunctionCheckBox_2      matlab.ui.control.CheckBox
        UIAxes4                       matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        c = 1
    end
    
    methods (Access = public)
        
        function addFlow(app)
            ft = {'Uniform', 'Source', 'Sink', 'Doublet', 'Vortex'};
            newFlow(1,1).FlowType = ft(randi([1 5]));
            newFlow(1,1).X = randi([-4 4]);
            newFlow(1,1).Y = randi([-2 2]);
            newFlow(1,1).Intensity = randi([1 10]);
            newFlow(1,1).Active = true;
            app.UITable.Data = [app.UITable.Data; struct2table(newFlow)];
        end
        
        function updatePlot(app)
            toll = 0.01;
            [x,y] = meshgrid(-5:toll:5,-2.5:toll:2.5);
            [x2,y2] = meshgrid(-5:5*toll:5,-2.5:5*toll:2.5);
            data = app.UITable.Data;
            %if app.JukowskitransformCheckBox.Value
            if 0
                [X,Y]=JukowskiTransform(app,x,y);
            else
                X = x;
                Y = y;
                X2 = x2;
                Y2 = y2;
            end
            Z = zeros(height(X), width(X));
            U = zeros(height(X2), width(X2));
            V = zeros(height(X2), width(X2));
            for c = 1:height(data)
                type = cellstr(data.FlowType(c));
                posX = data.X(c);
                posY = data.Y(c);
                Intensity = data.Intensity(c);
                switch type{:}
                    case 'Uniform'
                        app.UITable.Data.X(c) = nan;
                        app.UITable.Data.Y(c) = nan;
                        z = UniformFlow(app, X, Y, Intensity);
                        [u, v] = UniformFlowVel(app, X2, Y2, Intensity);
                    case 'Source'
                        z = SourceFlow(app, X, Y, posX, posY, Intensity);
                        [u, v] = SourceFlowVel(app, X2, Y2, posX, posY, Intensity);
                    case 'Sink'
                        z = SinkFlow(app, X, Y, posX, posY, Intensity);
                        [u, v] = SinkFlowVel(app, X2, Y2, posX, posY, Intensity);
                    case 'Doublet'
                        z = DoubletFlow(app, X, Y, posX, posY, Intensity);
                        [u, v] = DoubletFlowVel(app, X2, Y2, posX, posY, Intensity);
                    case 'Vortex'
                        z = VortexFlow(app, X, Y, posX, posY, Intensity);
                        [u, v] = VortexFlowVel(app, X2, Y2, posX, posY, Intensity);
                end
                if data.Active(c) == true
                    Z = Z + z;
                    U = U + u;
                    V = V + v;
                end
            end
            hold(app.UIAxes,'off')
            cla(app.UIAxes);
            if app.PressurefieldCheckBox.Value
                P = 0.5*1.225.*(U.^2+V.^2);
                c = colorbar(app.UIAxes);
                [~,h] = contourf(app.UIAxes, x2, y2, P, linspace(0, 100, 501));
                caxis(app.UIAxes,[0 100]);
                set(h,'LineColor','none');
                hold(app.UIAxes,'on');
            end
            if app.StreamfunctionCheckBox.Value
                contour(app.UIAxes, x, y, Z, linspace(-50,50,51), 'k');
                hold(app.UIAxes,'on')
            end
            if app.StreamplotCheckBox.Value
                streamslice(app.UIAxes, x2, y2, U, V, 2);
            end
        end
        
        function z = UniformFlow(~,~,y,intensity)
            z = intensity*y;
        end
        function [u,v] = UniformFlowVel(~,~,~,intensity)
            u = intensity;
            v = 0;
        end
        
        function z = SourceFlow(~,x,y,posX,posY,intensity)
            z = intensity*atan2(y-posY, x-posX)/(2*pi);
        end
        function [u,v] = SourceFlowVel(~,x,y,posX,posY,intensity)
            u = intensity/(2*pi)*(x-posX)./((x-posX).^2+(y-posY).^2);
            v = intensity/(2*pi)*(y-posY)./((x-posX).^2+(y-posY).^2);
        end
        
        function z = SinkFlow(~,x,y,posX,posY,intensity)
            z = -intensity*atan2(y-posY, x-posX)/(2*pi);
        end
        function [u,v] = SinkFlowVel(~,x,y,posX,posY,intensity)
            u = -intensity/(2*pi)*(x-posX)./((x-posX).^2+(y-posY).^2);
            v = -intensity/(2*pi)*(y-posY)./((x-posX).^2+(y-posY).^2);
        end
        
        function z = DoubletFlow(~,x,y,posX,posY,intensity)
            z = -intensity*(y-posY)./(2*pi*((x-posX).^2+(y-posY).^2));
        end
        function [u,v] = DoubletFlowVel(~,x,y,posX,posY,intensity)
            u = -intensity*((x-posX).^2-(y-posY).^2)./(2*pi*((x-posX).^2+(y-posY).^2).^2);
            v = -intensity*((x-posX).*(y-posY))./(2*pi*((x-posX).^2+(y-posY).^2).^2);
        end
        function z = VortexFlow(~,x,y,posX,posY,intensity)
            z = intensity/(4*pi)*log((x-posX).^2+(y-posY).^2);
        end
        function [u,v] = VortexFlowVel(~,x,y,posX,posY,intensity)
            u = intensity*(y-posY)./(2*pi*((x-posX).^2+(y-posY).^2));
            v = -intensity*(x-posX)./(2*pi*((x-posX).^2+(y-posY).^2));
        end
        
        
        function [xi, yi] = JukowskiTransform(~,x, y)
            z = x+1i.*y;
            epsilon = z;
            xi = real(epsilon);
            yi = imag(epsilon);
        end
    end
    
    methods (Access = private)
        
        function results = grad2rad(~, alpha)
            results = alpha * pi / 180;
        end
        
        function [J,f, z, z_circle, z_airfoil, L_str] = joukowski(app, v_inf,...
                alpha, x_cir, y_cir, radio, toll)
            
            r = radio;
            v = v_inf/v_inf;
            theta = grad2rad(app,alpha);
            s = x_cir + y_cir * 1i;
            
            % FLUID PARAMETER
            rho = 1.225;
            % TRANSFORMATION PARAMETER
            lambda = r-s;
            % CIRCULATION
            beta = (theta);
            k = 2*r*v*sin(beta);
            Gamma = k/(2*pi); %CIRCULATION

            %COMPLEX ASYMPTOTIC SPEED 
            w = v * exp(1i*theta);
            x = meshgrid(-3:.01:3);
            y = x';

            % COMPLEX PLANE
            z = x + 1i*y;

            for a = 1:length(x)
                for b = 1:length(y)
                    if abs(z(a,b)-s) <=  r - toll
                        z(a,b) = NaN;
                    end
                end
            end
            % AERODYNAMIC POTENTIAL
            f = w*(z) + (v*exp(-1i*theta)*r^2)./(z-s) + 1i*k*log(z);

            % JOUKOWSKI TRANSFORMATION, 
            J = z+lambda^2./z;

            %GRAPHIC - Circle and Joukowski Airfoil
            angle = 0:.1:2*pi;
            z_circle = r*(cos(angle)+1i*sin(angle)) + s;
            z_airfoil = z_circle+lambda^2./z_circle;
            
            % KUTTA JOUKOWSKI THEOREM
            L = v_inf*rho*Gamma;
            L_str = num2str(L);
            
        end
        
        function updateJoukowski(app)
            cla(app.UIAxes2)
            cla(app.UIAxes3)
            
            v_inf = app.VelocitySlider.Value;
            alpha = app.AoASlider.Value;
            x_cir = app.XEditField.Value;
            y_cir = app.YEditField.Value;
            radio = app.REditField.Value;
            toll = 10^app.TollSlider.Value;
            [J,f, z, z_circle, z_airfoil, L_str] = joukowski(app, v_inf,...
                alpha, x_cir, y_cir, radio, toll);
            
            ma = -5:.2:5;
            
            if app.PressurefieldCheckBox_2.Value
                [u, v] = velocityJoukowski(app, f);
                P = 0.5*1.225.*(u.^2+v.^2);
                [~, h] = contourf(app.UIAxes2, real(z), imag(z), P, linspace(0, 10*app.VelocitySlider.Value, 501));
                set(h,'LineColor','none');
                hold(app.UIAxes2, 'on')
            else
                contour(app.UIAxes2, real(z), imag(z), imag(f), ma, 'k')
            end
            hold(app.UIAxes2, 'on')
            fill(app.UIAxes2, real(z_circle), imag(z_circle),'y')
            set(app.UIAxes2,'DataAspectRatio',[1 1 1])
            axis(app.UIAxes2,[-3 3 -3 3])
            hold(app.UIAxes2, 'off')

            % title(strcat('Flow Around a Circle.   Lift:  ',L_str,'  [N/m]'));
            
            if app.PressurefieldCheckBox_2.Value
                [u, v] = velocityJoukowski(app, f);
                P = 0.5*1.225.*(u.^2+v.^2);
                [~, h] = contourf(app.UIAxes3, real(J), imag(J), P, linspace(0, 10*app.VelocitySlider.Value, 501));
                set(h,'LineColor','none');
                hold(app.UIAxes3, 'on')
            else
                contour(app.UIAxes3,real(J),imag(J),imag(f), ma, 'k')
            end
            hold(app.UIAxes3, 'on')
            fill(app.UIAxes3,real(z_airfoil),imag(z_airfoil),'y')
            set(app.UIAxes3,'DataAspectRatio',[1 1 1])
            axis(app.UIAxes3, [-3 3 -3 3])
            hold(app.UIAxes3, 'off')
            app.LLabel.Text = L_str;
            app.DLabel.Text = '0';
            % title(strcat('Flow Around the Corresponding Airfoil.   Lift:  ',L_str,'  [N/m]'));
        end
        
        function [u, v] = velocityJoukowski(~, psi)
            [psix, psiy] = gradient(psi, .01, .01);
            u = imag(psiy);
            v = imag(-psix);
        end
        
        
        function [m, p, t] = nacaParameters(app)
            naca = app.NACAairfoilEditField.Value;
            a = str2double(naca(1));
            b = str2double(naca(2));
            cd = str2double(naca(3:4));
            m = 0.01 * a;
            p = 0.1 * b;
            t = 0.01 * cd;
        end
        
        function [yc, xu, xl, yu, yl] = f(~, x, m, p, t, c)
            y1 = m/(p^2).*(2*p.*x-x.^2);
            y2 = m/((1-p)^2).*((1-2*p)+2*p.*x-x.^2);
            if p~=0
                yc = (x<p).*y1+(x>=p).*y2;
            else
                yc = y2;
            end
            yt = 5*t.*(0.2969.*(x.^(0.5))-0.126.*x-0.3516.*(x.^2) ...
                +0.2843.*(x.^3)-0.1015.*(x.^4));
            y3 = 2*m/(p.^2).*(p-x);
            y4 = 2*m/((1-p)^2).*(p-x);
            if p~=0
                dyc = (x<p).*y3+(x>=p).*y4;
            else
                dyc = y4;
            end
            theta = atan(dyc);
            xu = x - yt.*sin(theta);
            xl = x + yt.*sin(theta);
            yu = yc + yt.*cos(theta);
            yl = yc - yt.*cos(theta);
        end
        
        
        function panelPlot(app, m, p, t)
            x = 0:0.01:app.c;
            [yc, xu, xl, yu, yl] = f(app, x, m, p, t, app.c);
            cla(app.UIAxes4);
            plot(app.UIAxes4, x-app.c/2, yc, 'b');
            hold(app.UIAxes4, 'on');
            plot(app.UIAxes4, xu-app.c/2, yu, 'r');
            plot(app.UIAxes4, xl-app.c/2, yl, 'r');
            airfoil = polyshape({[-5,-5,5,5],[fliplr(xl-app.c/2) xu-app.c/2]}, ...
                {[-2,2,2,-2],[fliplr(yl) yu]});
            tr = triangulation(airfoil);
            model = createpde;
            tnodes = tr.Points';
            telements = tr.ConnectivityList';
            [G, mesh] = geometryFromMesh(model, tnodes, telements);
            pdegplot(model);
            generateMesh(model);
            set(0,'DefaultFigureVisible','on');
            pdemesh(model);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            t = table;
            t.FlowType = categorical({'Uniform'},{'Uniform', 'Source', ...
                'Sink', 'Doublet', 'Vortex'});
            t.X = nan;
            t.Y = nan;
            t.Intensity = 10;
            t.Active = true;
            app.UITable.Data = t;
            app.UITable.ColumnEditable = true;
            app.UIAxes.XLim = [-5 5];
            app.UIAxes.YLim = [-2.5 2.5];
            app.UIAxes4.XLim = [-5 5];
            app.UIAxes4.YLim = [-2 2];
            set(app.UIAxes4,'DataAspectRatio',[1 1 1]);
            set(app.UIAxes2,'DataAspectRatio',[1 1 1]);
            set(app.UIAxes3,'DataAspectRatio',[1 1 1]);
            app.StreamfunctionCheckBox.Value = true;
            updatePlot(app);
            app.REditField.Value = 1;
            app.XEditField.Value = 0.1;
            app.YEditField.Value = 0.1;
            updateJoukowski(app);
        end

        % Button pushed function: AddnewflowButton
        function AddnewflowButtonPushed(app, event)
            addFlow(app);
            updatePlot(app);
        end

        % Value changed function: StreamfunctionCheckBox
        function StreamfunctionCheckBoxValueChanged(app, event)
            updatePlot(app);
        end

        % Value changed function: StreamplotCheckBox
        function StreamplotCheckBoxValueChanged(app, event)
            updatePlot(app);
        end

        % Cell edit callback: UITable
        function UITableCellEdit(app, event)
            updatePlot(app);
        end

        % Button pushed function: ResetflowsButton
        function ResetflowsButtonPushed(app, event)
            startupFcn(app);
        end

        % Menu selected function: HalfRankineOvalMenu
        function HalfRankineOvalMenuSelected(app, event)
            t = table;
            t.FlowType = categorical({'Uniform'; 'Source'},{'Uniform', ...
                'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = [nan; -1];
            t.Y = [nan; 0];
            t.Intensity = [6; 20];
            t.Active = [true; true];
            app.UITable.Data = t;
            updatePlot(app);
        end

        % Menu selected function: RankineOvalMenu
        function RankineOvalMenuSelected(app, event)
            t = table;
            t.FlowType = categorical({'Uniform'; 'Source'; 'Sink'}, ...
                {'Uniform', 'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = [nan; -1; 1];
            t.Y = [nan; 0; 0];
            t.Intensity = [6; 20; 20];
            t.Active = [true; true; true];
            app.UITable.Data = t;
            updatePlot(app);
        end

        % Menu selected function: NonliftingflowoveracylinderMenu
        function NonliftingflowoveracylinderMenuSelected(app, event)
            t = table;
            t.FlowType = categorical({'Uniform'; 'Doublet'},{'Uniform', ...
                'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = [nan; 0];
            t.Y = [nan; 0];
            t.Intensity = [6; 20];
            t.Active = [true; true];
            app.UITable.Data = t;
            updatePlot(app);
            set(0,'DefaultFigureVisible','off');
            R = sqrt(20/(2*pi*6));
            plot(app.UIAxes, nsidedpoly(1000, 'Center', [0 0], ...
                 'Radius', R), 'FaceColor', 'r')
            alpha(0)
            hold off
        end

        % Menu selected function: LiftingflowoveracylinderMenu
        function LiftingflowoveracylinderMenuSelected(app, event)
            t = table;
            t.FlowType = categorical({'Uniform'; 'Doublet'; 'Vortex'}, ...
                {'Uniform', 'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = [nan; 0; 0];
            t.Y = [nan; 0; 0];
            t.Intensity = [6; 20; 30];
            t.Active = [true; true; true];
            app.UITable.Data = t;
            updatePlot(app);
            set(0,'DefaultFigureVisible','off');
            R = sqrt(20/(2*pi*6));
            plot(app.UIAxes, nsidedpoly(1000, 'Center', [0 0], ...
                'Radius', R), 'FaceColor', 'r')
            alpha(0)
            hold off
        end

        % Menu selected function: BathtubevortexMenu
        function BathtubevortexMenuSelected(app, event)
            t = table;
            t.FlowType = categorical({'Vortex'; 'Sink'},{'Uniform', ...
                'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = [0; 0];
            t.Y = [0; 0];
            t.Intensity = [30; 30];
            t.Active = [true; true];
            app.UITable.Data = t;
            updatePlot(app);
        end

        % Menu selected function: GroundeffectMenu
        function GroundeffectMenuSelected(app, event)
            t = table;
            t.FlowType = categorical({'Uniform'; 'Vortex'; 'Vortex'}, ...
                {'Uniform', 'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = [nan; 0; 0];
            t.Y = [nan; 1; -1];
            t.Intensity = [6; 40; -40];
            t.Active = [true; true; true];
            app.UITable.Data = t;
            updatePlot(app);
            X = [-5 5];
            Y = [0 0];
            plot(app.UIAxes, X, Y, "LineWidth", 2,"Color",'#804000');
        end

        % Menu selected function: FiniterowofvorticesMenu
        function FiniterowofvorticesMenuSelected(app, event)
            app.Panel.Visible = 'on';
        end

        % Button pushed function: AcceptButton
        function AcceptButtonPushed(app, event)
            app.Panel.Visible = 'off';
            t = table;
            ft = {};
            num = app.VortexnumberSpinner.Value;
            sep = app.SeparationEditField.Value;
            width = (num-1)*sep/2;
            for c = 1:num
                ft{end+1}='Vortex';
            end
            t.FlowType = categorical(transpose(ft), ...
                {'Uniform', 'Source', 'Sink', 'Doublet', 'Vortex'});
            t.X = transpose(linspace(-width,width,num));
            t.Y = zeros(num,1);
            t.Intensity = 25*app.IntensitySlider.Value/num*ones(num,1);
            t.Active = true(num,1);
            app.UITable.Data = t;
            updatePlot(app);
            set(0,'DefaultFigureVisible','off');
            X = [-width width];
            Y = [0 0];
            plot(app.UIAxes, X, Y, "LineWidth", 2,"Color",'r');
            alpha(.5)
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(app, event)
            app.Panel.Visible = 'off';
        end

        % Callback function
        function ButtonPushed(app, event)
            updateJoukowski(app);
        end

        % Value changed function: XEditField
        function XEditFieldValueChanged(app, event)
            updateJoukowski(app);
        end

        % Value changed function: YEditField
        function YEditFieldValueChanged(app, event)
            updateJoukowski(app);
        end

        % Value changed function: REditField
        function REditFieldValueChanged(app, event)
            updateJoukowski(app);
        end

        % Value changed function: AoASlider
        function AoASliderValueChanged(app, event)
            updateJoukowski(app);
        end

        % Value changed function: VelocitySlider
        function VelocitySliderValueChanged(app, event)
            updateJoukowski(app);
        end

        % Value changed function: TollSlider
        function TollSliderValueChanged(app, event)
            updateJoukowski(app);
        end

        % Value changed function: PressurefieldCheckBox
        function PressurefieldCheckBoxValueChanged(app, event)
            updatePlot(app);
        end

        % Value changed function: NACAairfoilEditField
        function NACAairfoilEditFieldValueChanged(app, event)
            [m, p, t] = nacaParameters(app);
            panelPlot(app, m, p, t);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create Menu
            app.Menu = uimenu(app.UIFigure);
            app.Menu.Text = 'Menu';

            % Create BasicconfigurationsMenu
            app.BasicconfigurationsMenu = uimenu(app.UIFigure);
            app.BasicconfigurationsMenu.Text = 'Basic configurations     ';

            % Create HalfRankineOvalMenu
            app.HalfRankineOvalMenu = uimenu(app.BasicconfigurationsMenu);
            app.HalfRankineOvalMenu.MenuSelectedFcn = createCallbackFcn(app, @HalfRankineOvalMenuSelected, true);
            app.HalfRankineOvalMenu.Text = 'Half Rankine Oval';

            % Create RankineOvalMenu
            app.RankineOvalMenu = uimenu(app.BasicconfigurationsMenu);
            app.RankineOvalMenu.MenuSelectedFcn = createCallbackFcn(app, @RankineOvalMenuSelected, true);
            app.RankineOvalMenu.Text = 'Rankine Oval';

            % Create NonliftingflowoveracylinderMenu
            app.NonliftingflowoveracylinderMenu = uimenu(app.BasicconfigurationsMenu);
            app.NonliftingflowoveracylinderMenu.MenuSelectedFcn = createCallbackFcn(app, @NonliftingflowoveracylinderMenuSelected, true);
            app.NonliftingflowoveracylinderMenu.Text = 'Non-lifting flow over a cylinder';

            % Create LiftingflowoveracylinderMenu
            app.LiftingflowoveracylinderMenu = uimenu(app.BasicconfigurationsMenu);
            app.LiftingflowoveracylinderMenu.MenuSelectedFcn = createCallbackFcn(app, @LiftingflowoveracylinderMenuSelected, true);
            app.LiftingflowoveracylinderMenu.Text = 'Lifting flow over a cylinder';

            % Create BathtubevortexMenu
            app.BathtubevortexMenu = uimenu(app.BasicconfigurationsMenu);
            app.BathtubevortexMenu.MenuSelectedFcn = createCallbackFcn(app, @BathtubevortexMenuSelected, true);
            app.BathtubevortexMenu.Text = 'Bath-tube vortex';

            % Create GroundeffectMenu
            app.GroundeffectMenu = uimenu(app.BasicconfigurationsMenu);
            app.GroundeffectMenu.MenuSelectedFcn = createCallbackFcn(app, @GroundeffectMenuSelected, true);
            app.GroundeffectMenu.Text = 'Ground effect';

            % Create FiniterowofvorticesMenu
            app.FiniterowofvorticesMenu = uimenu(app.BasicconfigurationsMenu);
            app.FiniterowofvorticesMenu.MenuSelectedFcn = createCallbackFcn(app, @FiniterowofvorticesMenuSelected, true);
            app.FiniterowofvorticesMenu.Text = 'Finite row of vortices';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [2 1 640 480];

            % Create IntroTab
            app.IntroTab = uitab(app.TabGroup);
            app.IntroTab.Title = 'Intro';
            app.IntroTab.BackgroundColor = [1 1 1];
            app.IntroTab.Scrollable = 'on';

            % Create Image5
            app.Image5 = uiimage(app.IntroTab);
            app.Image5.Position = [149 204 344 204];
            app.Image5.ImageSource = 'aerodinamica_aviones.png';

            % Create Image6
            app.Image6 = uiimage(app.IntroTab);
            app.Image6.Position = [217 150 206 100];
            app.Image6.ImageSource = 'escudo.svg';

            % Create AERODYNAMICFLOWCONFIGURATIONLabel
            app.AERODYNAMICFLOWCONFIGURATIONLabel = uilabel(app.IntroTab);
            app.AERODYNAMICFLOWCONFIGURATIONLabel.FontName = 'Times New Roman';
            app.AERODYNAMICFLOWCONFIGURATIONLabel.FontSize = 26;
            app.AERODYNAMICFLOWCONFIGURATIONLabel.FontWeight = 'bold';
            app.AERODYNAMICFLOWCONFIGURATIONLabel.FontAngle = 'italic';
            app.AERODYNAMICFLOWCONFIGURATIONLabel.Position = [64 406 512 32];
            app.AERODYNAMICFLOWCONFIGURATIONLabel.Text = 'AERODYNAMIC FLOW CONFIGURATION';

            % Create Authors
            app.Authors = uilabel(app.IntroTab);
            app.Authors.HorizontalAlignment = 'center';
            app.Authors.FontWeight = 'bold';
            app.Authors.Position = [251 89 139 42];
            app.Authors.Text = {'Authors:'; 'Mario Álvarez Redondo'; 'Gabriel Valles Valverde'};

            % Create Date05102021Label
            app.Date05102021Label = uilabel(app.IntroTab);
            app.Date05102021Label.HorizontalAlignment = 'center';
            app.Date05102021Label.Position = [271 58 98 22];
            app.Date05102021Label.Text = 'Date: 05/10/2021';

            % Create BasicflowanalysisTab
            app.BasicflowanalysisTab = uitab(app.TabGroup);
            app.BasicflowanalysisTab.Title = 'Basic flow analysis';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.BasicflowanalysisTab);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {22, '2x', '1x', 22};
            app.GridLayout.RowSpacing = 12.1666666666667;
            app.GridLayout.Padding = [10 12.1666666666667 10 12.1666666666667];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            app.UIAxes.Toolbar.Visible = 'off';
            app.UIAxes.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
            app.UIAxes.XTickLabel = {'-5'; '-4'; '-3'; '-2'; '-1'; '0'; '1'; '2'; '3'; '4'; '5'};
            app.UIAxes.YTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
            app.UIAxes.YTickLabel = {'-5'; '-4'; '-3'; '-2'; '-1'; '0'; '1'; '2'; '3'; '4'; '5'};
            app.UIAxes.BoxStyle = 'full';
            app.UIAxes.ClippingStyle = 'rectangle';
            app.UIAxes.Box = 'on';
            app.UIAxes.Layout.Row = 2;
            app.UIAxes.Layout.Column = [2 5];

            % Create StreamfunctionCheckBox
            app.StreamfunctionCheckBox = uicheckbox(app.GridLayout);
            app.StreamfunctionCheckBox.ValueChangedFcn = createCallbackFcn(app, @StreamfunctionCheckBoxValueChanged, true);
            app.StreamfunctionCheckBox.Text = {'Stream function'; ''};
            app.StreamfunctionCheckBox.Layout.Row = 1;
            app.StreamfunctionCheckBox.Layout.Column = [5 6];

            % Create StreamplotCheckBox
            app.StreamplotCheckBox = uicheckbox(app.GridLayout);
            app.StreamplotCheckBox.ValueChangedFcn = createCallbackFcn(app, @StreamplotCheckBoxValueChanged, true);
            app.StreamplotCheckBox.Text = 'Stream plot';
            app.StreamplotCheckBox.Layout.Row = 1;
            app.StreamplotCheckBox.Layout.Column = 4;

            % Create ResetflowsButton
            app.ResetflowsButton = uibutton(app.GridLayout, 'push');
            app.ResetflowsButton.ButtonPushedFcn = createCallbackFcn(app, @ResetflowsButtonPushed, true);
            app.ResetflowsButton.Layout.Row = 4;
            app.ResetflowsButton.Layout.Column = 5;
            app.ResetflowsButton.Text = 'Reset flows';

            % Create AddnewflowButton
            app.AddnewflowButton = uibutton(app.GridLayout, 'push');
            app.AddnewflowButton.ButtonPushedFcn = createCallbackFcn(app, @AddnewflowButtonPushed, true);
            app.AddnewflowButton.Layout.Row = 4;
            app.AddnewflowButton.Layout.Column = 4;
            app.AddnewflowButton.Text = 'Add new flow';

            % Create UITable
            app.UITable = uitable(app.GridLayout);
            app.UITable.ColumnName = {'Flow Type'; 'X'; 'Y'; 'Intensity'; 'Active'};
            app.UITable.RowName = {};
            app.UITable.CellEditCallback = createCallbackFcn(app, @UITableCellEdit, true);
            app.UITable.Layout.Row = 3;
            app.UITable.Layout.Column = [2 5];

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Title = 'Panel';
            app.Panel.Visible = 'off';
            app.Panel.Layout.Row = 2;
            app.Panel.Layout.Column = [3 4];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel);
            app.GridLayout2.ColumnWidth = {54, 30, 55, '1x'};
            app.GridLayout2.RowHeight = {22, 22, 22, '1x', 22};
            app.GridLayout2.ColumnSpacing = 6.2;
            app.GridLayout2.Padding = [6.2 10 6.2 10];

            % Create AcceptButton
            app.AcceptButton = uibutton(app.GridLayout2, 'push');
            app.AcceptButton.ButtonPushedFcn = createCallbackFcn(app, @AcceptButtonPushed, true);
            app.AcceptButton.Layout.Row = 5;
            app.AcceptButton.Layout.Column = [1 2];
            app.AcceptButton.Text = 'Accept';

            % Create SeparationEditFieldLabel
            app.SeparationEditFieldLabel = uilabel(app.GridLayout2);
            app.SeparationEditFieldLabel.HorizontalAlignment = 'right';
            app.SeparationEditFieldLabel.Layout.Row = 2;
            app.SeparationEditFieldLabel.Layout.Column = [1 2];
            app.SeparationEditFieldLabel.Text = 'Separation';

            % Create SeparationEditField
            app.SeparationEditField = uieditfield(app.GridLayout2, 'numeric');
            app.SeparationEditField.Limits = [0 1];
            app.SeparationEditField.Layout.Row = 2;
            app.SeparationEditField.Layout.Column = [3 4];
            app.SeparationEditField.Value = 0.01;

            % Create CloseButton
            app.CloseButton = uibutton(app.GridLayout2, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Layout.Row = 5;
            app.CloseButton.Layout.Column = [3 4];
            app.CloseButton.Text = 'Close';

            % Create VortexnumberSpinnerLabel
            app.VortexnumberSpinnerLabel = uilabel(app.GridLayout2);
            app.VortexnumberSpinnerLabel.HorizontalAlignment = 'right';
            app.VortexnumberSpinnerLabel.Layout.Row = 1;
            app.VortexnumberSpinnerLabel.Layout.Column = [1 2];
            app.VortexnumberSpinnerLabel.Text = 'Vortex number';

            % Create VortexnumberSpinner
            app.VortexnumberSpinner = uispinner(app.GridLayout2);
            app.VortexnumberSpinner.Limits = [1 500];
            app.VortexnumberSpinner.Layout.Row = 1;
            app.VortexnumberSpinner.Layout.Column = [3 4];
            app.VortexnumberSpinner.Value = 200;

            % Create IntensitySlider
            app.IntensitySlider = uislider(app.GridLayout2);
            app.IntensitySlider.Limits = [1 8];
            app.IntensitySlider.MajorTicks = [1 8];
            app.IntensitySlider.MajorTickLabels = {'Min ', 'Max'};
            app.IntensitySlider.Layout.Row = 3;
            app.IntensitySlider.Layout.Column = [3 4];
            app.IntensitySlider.Value = 1;

            % Create IntensitySliderLabel
            app.IntensitySliderLabel = uilabel(app.GridLayout2);
            app.IntensitySliderLabel.HorizontalAlignment = 'right';
            app.IntensitySliderLabel.Layout.Row = 3;
            app.IntensitySliderLabel.Layout.Column = [1 2];
            app.IntensitySliderLabel.Text = 'Intensity';

            % Create PressurefieldCheckBox
            app.PressurefieldCheckBox = uicheckbox(app.GridLayout);
            app.PressurefieldCheckBox.ValueChangedFcn = createCallbackFcn(app, @PressurefieldCheckBoxValueChanged, true);
            app.PressurefieldCheckBox.Text = 'Pressure field';
            app.PressurefieldCheckBox.Layout.Row = 1;
            app.PressurefieldCheckBox.Layout.Column = 3;

            % Create JoukowskitransformationTab
            app.JoukowskitransformationTab = uitab(app.TabGroup);
            app.JoukowskitransformationTab.Title = 'Joukowski transformation';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.JoukowskitransformationTab);
            app.GridLayout3.ColumnWidth = {'1x'};

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.GridLayout3);
            app.GridLayout4.ColumnWidth = {'1x', '2x'};
            app.GridLayout4.RowHeight = {'1x'};
            app.GridLayout4.Layout.Row = 2;
            app.GridLayout4.Layout.Column = 1;

            % Create GridLayout4_2
            app.GridLayout4_2 = uigridlayout(app.GridLayout4);
            app.GridLayout4_2.ColumnWidth = {'1x', '3x', '0.5x'};
            app.GridLayout4_2.RowHeight = {22, 22, 22, 22, '1x'};
            app.GridLayout4_2.Layout.Row = 1;
            app.GridLayout4_2.Layout.Column = 1;

            % Create XEditFieldLabel
            app.XEditFieldLabel = uilabel(app.GridLayout4_2);
            app.XEditFieldLabel.HorizontalAlignment = 'right';
            app.XEditFieldLabel.Layout.Row = 1;
            app.XEditFieldLabel.Layout.Column = 1;
            app.XEditFieldLabel.Text = 'X';

            % Create XEditField
            app.XEditField = uieditfield(app.GridLayout4_2, 'numeric');
            app.XEditField.ValueChangedFcn = createCallbackFcn(app, @XEditFieldValueChanged, true);
            app.XEditField.Layout.Row = 1;
            app.XEditField.Layout.Column = 2;

            % Create YEditFieldLabel
            app.YEditFieldLabel = uilabel(app.GridLayout4_2);
            app.YEditFieldLabel.HorizontalAlignment = 'right';
            app.YEditFieldLabel.Layout.Row = 2;
            app.YEditFieldLabel.Layout.Column = 1;
            app.YEditFieldLabel.Text = 'Y';

            % Create YEditField
            app.YEditField = uieditfield(app.GridLayout4_2, 'numeric');
            app.YEditField.ValueChangedFcn = createCallbackFcn(app, @YEditFieldValueChanged, true);
            app.YEditField.Layout.Row = 2;
            app.YEditField.Layout.Column = 2;

            % Create REditFieldLabel
            app.REditFieldLabel = uilabel(app.GridLayout4_2);
            app.REditFieldLabel.HorizontalAlignment = 'right';
            app.REditFieldLabel.Layout.Row = 3;
            app.REditFieldLabel.Layout.Column = 1;
            app.REditFieldLabel.Text = 'R';

            % Create REditField
            app.REditField = uieditfield(app.GridLayout4_2, 'numeric');
            app.REditField.ValueChangedFcn = createCallbackFcn(app, @REditFieldValueChanged, true);
            app.REditField.Layout.Row = 3;
            app.REditField.Layout.Column = 2;

            % Create Image7
            app.Image7 = uiimage(app.GridLayout4_2);
            app.Image7.Layout.Row = 5;
            app.Image7.Layout.Column = [1 3];
            app.Image7.ImageSource = 'logo_ingenierias_retina.png';

            % Create PressurefieldCheckBox_2
            app.PressurefieldCheckBox_2 = uicheckbox(app.GridLayout4_2);
            app.PressurefieldCheckBox_2.Text = 'Pressure field';
            app.PressurefieldCheckBox_2.Layout.Row = 4;
            app.PressurefieldCheckBox_2.Layout.Column = 2;

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.GridLayout4);
            app.GridLayout5.ColumnWidth = {70, '1x'};
            app.GridLayout5.RowHeight = {'1x', '1x', '1x', '1.2x'};
            app.GridLayout5.Layout.Row = 1;
            app.GridLayout5.Layout.Column = 2;

            % Create AoASliderLabel
            app.AoASliderLabel = uilabel(app.GridLayout5);
            app.AoASliderLabel.Interpreter = 'latex';
            app.AoASliderLabel.HorizontalAlignment = 'right';
            app.AoASliderLabel.Layout.Row = 1;
            app.AoASliderLabel.Layout.Column = 1;
            app.AoASliderLabel.Text = '$\alpha$';

            % Create AoASlider
            app.AoASlider = uislider(app.GridLayout5);
            app.AoASlider.Limits = [-15 15];
            app.AoASlider.MajorTicks = [-15 -10 -5 0 5 10 15];
            app.AoASlider.ValueChangedFcn = createCallbackFcn(app, @AoASliderValueChanged, true);
            app.AoASlider.Layout.Row = 1;
            app.AoASlider.Layout.Column = 2;

            % Create VelocitySliderLabel
            app.VelocitySliderLabel = uilabel(app.GridLayout5);
            app.VelocitySliderLabel.Interpreter = 'latex';
            app.VelocitySliderLabel.HorizontalAlignment = 'right';
            app.VelocitySliderLabel.Layout.Row = 2;
            app.VelocitySliderLabel.Layout.Column = 1;
            app.VelocitySliderLabel.Text = '$V_{\infty}$';

            % Create VelocitySlider
            app.VelocitySlider = uislider(app.GridLayout5);
            app.VelocitySlider.Limits = [10 100];
            app.VelocitySlider.MajorTicks = [10 20 30 40 50 60 70 80 90 100];
            app.VelocitySlider.MajorTickLabels = {'10', '20', '30', '40', '50', '60', '70', '80', '90', '100', '', ''};
            app.VelocitySlider.ValueChangedFcn = createCallbackFcn(app, @VelocitySliderValueChanged, true);
            app.VelocitySlider.Layout.Row = 2;
            app.VelocitySlider.Layout.Column = 2;
            app.VelocitySlider.Value = 10;

            % Create TollSlider
            app.TollSlider = uislider(app.GridLayout5);
            app.TollSlider.Limits = [-5 -3];
            app.TollSlider.MajorTicks = [-5 -4 -3];
            app.TollSlider.MajorTickLabels = {'1e-5', '1e-4', '1e-3', ''};
            app.TollSlider.ValueChangedFcn = createCallbackFcn(app, @TollSliderValueChanged, true);
            app.TollSlider.Layout.Row = 3;
            app.TollSlider.Layout.Column = 2;
            app.TollSlider.Value = -4;

            % Create TollSliderLabel
            app.TollSliderLabel = uilabel(app.GridLayout5);
            app.TollSliderLabel.HorizontalAlignment = 'right';
            app.TollSliderLabel.Layout.Row = 3;
            app.TollSliderLabel.Layout.Column = 1;
            app.TollSliderLabel.Text = 'Toll';

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout5);
            app.GridLayout7.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout7.RowHeight = {'1x'};
            app.GridLayout7.Layout.Row = 4;
            app.GridLayout7.Layout.Column = [1 2];

            % Create LiftLabel
            app.LiftLabel = uilabel(app.GridLayout7);
            app.LiftLabel.Interpreter = 'latex';
            app.LiftLabel.HorizontalAlignment = 'right';
            app.LiftLabel.FontSize = 14;
            app.LiftLabel.FontWeight = 'bold';
            app.LiftLabel.Layout.Row = 1;
            app.LiftLabel.Layout.Column = 1;
            app.LiftLabel.Text = '$Lift$';

            % Create NmLabel
            app.NmLabel = uilabel(app.GridLayout7);
            app.NmLabel.Layout.Row = 1;
            app.NmLabel.Layout.Column = 3;
            app.NmLabel.Text = 'N/m';

            % Create DragLabel
            app.DragLabel = uilabel(app.GridLayout7);
            app.DragLabel.Interpreter = 'latex';
            app.DragLabel.HorizontalAlignment = 'right';
            app.DragLabel.FontSize = 14;
            app.DragLabel.FontWeight = 'bold';
            app.DragLabel.Layout.Row = 1;
            app.DragLabel.Layout.Column = 4;
            app.DragLabel.Text = '$Drag$';

            % Create NmLabel_2
            app.NmLabel_2 = uilabel(app.GridLayout7);
            app.NmLabel_2.Layout.Row = 1;
            app.NmLabel_2.Layout.Column = 6;
            app.NmLabel_2.Text = 'N/m';

            % Create LLabel
            app.LLabel = uilabel(app.GridLayout7);
            app.LLabel.HorizontalAlignment = 'right';
            app.LLabel.Layout.Row = 1;
            app.LLabel.Layout.Column = 2;
            app.LLabel.Text = 'L';

            % Create DLabel
            app.DLabel = uilabel(app.GridLayout7);
            app.DLabel.HorizontalAlignment = 'right';
            app.DLabel.Layout.Row = 1;
            app.DLabel.Layout.Column = 5;
            app.DLabel.Text = 'D';

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout3);
            app.GridLayout6.RowHeight = {'1x'};
            app.GridLayout6.Layout.Row = 1;
            app.GridLayout6.Layout.Column = 1;

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout6);
            title(app.UIAxes2, 'Title')
            app.UIAxes2.Layout.Row = 1;
            app.UIAxes2.Layout.Column = 1;

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.GridLayout6);
            title(app.UIAxes3, 'Title')
            app.UIAxes3.Layout.Row = 1;
            app.UIAxes3.Layout.Column = 2;

            % Create PanelmethodanalysisTab
            app.PanelmethodanalysisTab = uitab(app.TabGroup);
            app.PanelmethodanalysisTab.Title = 'Panel method analysis';

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.PanelmethodanalysisTab);
            app.GridLayout10.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout10.RowHeight = {22, '2x', 22, 22};
            app.GridLayout10.RowSpacing = 12.1666666666667;
            app.GridLayout10.Padding = [10 12.1666666666667 10 12.1666666666667];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.GridLayout10);
            app.UIAxes4.Toolbar.Visible = 'off';
            app.UIAxes4.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
            app.UIAxes4.XTickLabel = {'-5'; '-4'; '-3'; '-2'; '-1'; '0'; '1'; '2'; '3'; '4'; '5'};
            app.UIAxes4.YTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
            app.UIAxes4.YTickLabel = {'-5'; '-4'; '-3'; '-2'; '-1'; '0'; '1'; '2'; '3'; '4'; '5'};
            app.UIAxes4.BoxStyle = 'full';
            app.UIAxes4.ClippingStyle = 'rectangle';
            app.UIAxes4.Box = 'on';
            app.UIAxes4.Layout.Row = 2;
            app.UIAxes4.Layout.Column = [1 6];

            % Create StreamfunctionCheckBox_2
            app.StreamfunctionCheckBox_2 = uicheckbox(app.GridLayout10);
            app.StreamfunctionCheckBox_2.Text = {'Stream function'; ''};
            app.StreamfunctionCheckBox_2.Layout.Row = 1;
            app.StreamfunctionCheckBox_2.Layout.Column = [5 6];

            % Create StreamplotCheckBox_2
            app.StreamplotCheckBox_2 = uicheckbox(app.GridLayout10);
            app.StreamplotCheckBox_2.Text = 'Stream plot';
            app.StreamplotCheckBox_2.Layout.Row = 1;
            app.StreamplotCheckBox_2.Layout.Column = 4;

            % Create ComputepropertiesButton
            app.ComputepropertiesButton = uibutton(app.GridLayout10, 'push');
            app.ComputepropertiesButton.Layout.Row = 4;
            app.ComputepropertiesButton.Layout.Column = [3 4];
            app.ComputepropertiesButton.Text = 'Compute properties';

            % Create PressurefieldCheckBox_3
            app.PressurefieldCheckBox_3 = uicheckbox(app.GridLayout10);
            app.PressurefieldCheckBox_3.Text = 'Pressure field';
            app.PressurefieldCheckBox_3.Layout.Row = 1;
            app.PressurefieldCheckBox_3.Layout.Column = 3;

            % Create NACAairfoilEditFieldLabel
            app.NACAairfoilEditFieldLabel = uilabel(app.GridLayout10);
            app.NACAairfoilEditFieldLabel.HorizontalAlignment = 'right';
            app.NACAairfoilEditFieldLabel.Layout.Row = 3;
            app.NACAairfoilEditFieldLabel.Layout.Column = 3;
            app.NACAairfoilEditFieldLabel.Text = 'NACA airfoil';

            % Create NACAairfoilEditField
            app.NACAairfoilEditField = uieditfield(app.GridLayout10, 'text');
            app.NACAairfoilEditField.ValueChangedFcn = createCallbackFcn(app, @NACAairfoilEditFieldValueChanged, true);
            app.NACAairfoilEditField.Layout.Row = 3;
            app.NACAairfoilEditField.Layout.Column = 4;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
