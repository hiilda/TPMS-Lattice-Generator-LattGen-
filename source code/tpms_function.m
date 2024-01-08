function f = tpms_function(x1,y1,z1,lattice_type,custom_function)

X = 2*pi*x1;
Y = 2*pi*y1;
Z = 2*pi*z1;

switch lattice_type
    case "Primitive Schwartz Surface"
        f = cos(X)+cos(Y)+cos(Z);
    case "Double P Surface"
        f = 0.5*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))+0.2*(cos(2*X)+cos(2*Y)+cos(2*Z));
    case "Gyroid Surface"
        f = sin(X).*cos(Y) + sin(Z).*cos(X) + sin(Y).*cos(Z);
    case "Double Gyroid Surface"
        f = 2.75.*(sin(2*X).*sin(Z).*cos(Y) + sin(2*Y).*sin(X).*cos(Z) + sin(2*Z).*sin(Y).*cos(X))-(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X));
    case "Diamond Surface"
        f = (cos(X).*cos(Y).*cos(Z)) - (sin(X).*sin(Y).*sin(Z));
    case "Double Diamond 1 Surface"
        f = sin(2*X).*sin(2*Y)+sin(2*Y).*sin(2*Z)+sin(2*X).*sin(2*Z)+cos(2*X).*cos(2*Y).*cos(2*Z);
    case "Double Diamond 2 Surface"
        f = cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*X).*cos(2*Z)+sin(2*X).*sin(2*Y).*sin(2*Z);
    case "IWP Surface"
        f = 2*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-(cos(2*X)+cos(2*Y)+cos(2*Z));
    case "FRD Prime Surface"
        f = 4*(cos(X).*cos(Y).*cos(Z))-(cos(2*X).*cos(2*Y)+cos(2*X).*cos(2*Z)+cos(2*Y).*cos(2*Z));
    case "I2-Y Surface"
        f = -2*(sin(2*X).*cos(Y).*sin(Z) + sin(X).*sin(2*Y).*cos(Z) + cos(X).*sin(Y).*sin(2*Z)) + cos(2*X).*cos(2*Y) + cos(2*Y).*cos(2*Z) + cos(2*X).*cos(2*Z);
    case "G Prime 1 Surface"
        f = 5*(sin(2*X).*sin(Z).*cos(Y)+sin(2*Y).*sin(X).*cos(Z)+sin(2*Z).*sin(Y).*cos(X))+(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X));
    case "G Prime 2 Surface"
        f = (sin(2*X).*sin(Z).*cos(Y)+sin(2*Y).*sin(X).*cos(Z)+sin(2*Z).*sin(Y).*cos(X))+0.32;
    case "Neovius Surface"
        f = 3*(cos(X)+cos(Y)+cos(Z))+4*(cos(X).*cos(Y).*cos(Z));
    case "Lidinoid Surface"
        f = (sin(2*X).*sin(Z).*cos(Y)+sin(2*Y).*sin(X).*cos(Z)+sin(2*Z).*sin(Y).*cos(X))-(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))+0.3;
    case "D Prime Surface"
        f = 0.5*(sin(X).*sin(Y).*sin(Z)+cos(X).*cos(Y).*cos(Z))-0.5*(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))-0.2;
    case "KP Surface"
        f = 0.3*(cos(X)+cos(Y)+cos(Z))+0.3*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-0.4*(cos(2*X)+cos(2*Y)+cos(2*Z))+0.2;
    case "S Surface"
        f = cos(2*X).*sin(Y).*cos(Z)+cos(2*Y).*sin(Z).*cos(X)+cos(2*Z).*sin(X).*cos(Y)-0.4;
    case "Fischer-Koch Surface"
        f = (cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-(cos(2*X)+cos(2*Y)+cos(2*Z));
    case "Bionic Bone 1 Surface"
        f = 20*(cos(X).*sin(Y)+cos(Y).*sin(Z)+cos(Z).*sin(X))-0.5*(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))-4;
    case "Bionic Bone 2 Surface"
        f = 10*(cos(X).*sin(Y)+cos(Y).*sin(Z)+cos(Z).*sin(X))-2*(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))-12;
    case "Octo 1 Surface"
        f = 4*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-2.8*(cos(X).*cos(Y).*cos(Z))+(cos(X)+cos(Y)+cos(Z))+1.5;
    case "Octo 2 Surface"
        f = 4*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X))-3*(cos(X)+cos(Y)+cos(Z))+2.4;
    case "PN Surface"
        f = 0.6*(cos(X).*cos(Y).*cos(Z))+0.4*(cos(X)+cos(Y)+cos(Z))+0.2*(cos(2*X).*cos(2*Y).*cos(2*Z))+0.2*(cos(2*X)+cos(2*Y)+cos(2*Z))+0.1*(cos(3*X)+cos(3*Y)+cos(3*Z))+0.2*(cos(X).*cos(Y)+cos(Y).*cos(Z)+cos(Z).*cos(X));
    case "C(Y) Surface"
        f = sin(X).*sin(Y).*sin(Z)+sin(2*X).*sin(Y)+sin(2*Y).*sin(Z)+sin(2*Z).*sin(X)-cos(X).*cos(Y).*cos(Z)+sin(2*X).*cos(Z)+sin(2*Y).*cos(X)+sin(2*Z).*cos(Y);
    case "Black D Surface"
        f = sin(X).*sin(Y).*sin(Z)+sin(X).*cos(Y).*cos(Z)+cos(X).*sin(Y).*cos(Z)+cos(X).*cos(Y).*sin(Z);
    case "FRD Surface"
        f = 8*(cos(X).*cos(Y).*cos(Z))+cos(2*X).*cos(2*Y).*cos(2*Z)-cos(2*X).*cos(2*Y)-cos(2*Y).*cos(2*Z)-cos(2*Z).*cos(2*X);
    case "Diamond D Surface"
        f = sin(X).*sin(Y).*sin(Z)+sin(X).*cos(Y).*cos(Z)+cos(X).*sin(Y).*cos(Z)+cos(X).*cos(Y).*sin(Z);
    case "D Surface"
        f = cos(X).*cos(Y).*cos(Z)-sin(X).*sin(Y).*sin(Z);
    case "Split P Surface"
        f = 1.1*(sin(2*X).*sin(Z).*cos(Y)+sin(2*Y).*sin(X).*cos(Z)+sin(2*Z).*sin(Y).*cos(X))-0.2*(cos(2*X).*cos(2*Y)+cos(2*Y).*cos(2*Z)+cos(2*Z).*cos(2*X))-0.4*(cos(X)+cos(Y)+cos(Z));
    case "Custom Function"
        if ischar(custom_function) == 1
            f = eval(custom_function);
        else
            f = custom_function;
        end
end