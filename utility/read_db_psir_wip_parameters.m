function out = read_db_psir_wip_parameters(filename);

[patientID, protocol] = read_ismrmrd_protocol(filename);

seqname = protocol.tSequenceFileName;

% if ~isempty(findstr(seqname,'PSIR_T2'))

    out.seqname = seqname;
    out.T1_myo = protocol.sWipMemBlock.alFree_{28};
    out.T1_blood = protocol.sWipMemBlock.alFree_{29};
    out.TD1 = protocol.sWipMemBlock.alFree_{30};
    out.TD2 = protocol.sWipMemBlock.alFree_{31};
    out.delta = protocol.sWipMemBlock.alFree_{32};
    out.TE = protocol.sWipMemBlock.alFree_{26};
    
% else
%     out = [];
% end
