classdef dipolePairAndMeasure
    properties
        from = pr.dipole;
        to = pr.dipole;
        toDipoleId % since dipoles are shared across pairs, this uniquely identified shared dipoles.
        fromDipoleId       
        sessionId % to calculate how many session contribute to a certain graph edge       
        linearizedMeasure = []; % N x M matrix, N is the given measure (ERP, ERSP...) in a linearized manner, M is the number of items (usually dipoles).
        measureLabel = ''; % the text label describing the measure, for example 'ERSP' or 'ERP'. This is set in the child classes;
    end;
        
    methods
    end;
end