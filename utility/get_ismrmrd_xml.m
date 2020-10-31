function hdr = get_ismrmrd_xml(h5name)
% hdr = get_ismrmrd_xml(h5name)

dset = ismrmrd.Dataset(h5name, 'dataset');
hdr = ismrmrd.xml.deserialize(dset.readxml());
