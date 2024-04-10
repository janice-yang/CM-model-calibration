%{
AbsCM is an abstract class representing an abstract cardiomyocyte. It inherits
from the superclass handle. It has concrete subclasses (eg. paci18, ohara,
kernik19 etc) that represent the different types of cells that could be
instantiated. AbsCM lays out properties of all CMs and methods that ought to be 
common to all CMs.
%}

classdef (Abstract) AbsCM < handle % We use handle because we want to treat
    % CMs as objects rather than values. Changes to that CM should be
    % reflected in all pointers to the CM.
    properties (Abstract)
        name % a string representation of the cell
        ODEModel % the 'dydt' function for the CM (e.g. @f_Paci2018, 
        % @f_Kernik2019, @f_ORd, etc)
        geometry % contains the geometric parameters of the CM (e.g. V_SR, Vc, 
        % Cm, jSR, etc.)
        state % contains the initial state variables as well as how they change 
        % alongside a vector of time points. Also contains the currents as they
        % change for alongside these time points.
        parameters % contains the baseline parameters vector as well as a vector 
        % of scaling factors if required
        conductances % contains the baseline conductances vector as well as a 
        % vector of scaling factors if required and details of the drug
        % application protocol if required
        protocol % contains the pacing protocol
    end
    methods (Abstract)
        setUpPacingProtocol(CM) % function that will create the pacing protocol
        % for the CM and store it in the CM's protocol struct
        
        odeSolver(CM) % function that will take the CM and run the ODE solver 
        % according to the CM's protocol using the CM's parameters 
        
        getCurrents(CM) % function that will produce & store a matrix 
        % of each of the currents & E_ion values over time. 
        
        %makePlots(CM) % function that will produce plots associated with the CM
    end

end