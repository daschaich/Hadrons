import h5py


print('meson correlator - ScalarSink')
f = h5py.File('mesons/pt_ll.1500.h5','r')
for ct in f['meson']['meson_0']['corr']:
    print(ct)
print()


print('scMat meson correlator - ScalarSink')
f = h5py.File('mesons/pion.1500.h5','r')
for ct in f['pion']['corr']:
    print(ct)
print()


#print('meson correlator - sources')
#f = h5py.File('mesons/pt_sourcell.1500.h5','r')
#for ct in f['meson']['meson_0']['corr']:
#    print(ct)
#print()

#print('meson correlator - Slice Propagator')
#f = h5py.File('mesons/pt_smeared_ll.1500.h5','r')
#for ct in f['meson']['meson_0']['corr']:
#    print(ct)
#print()


print('s0 in pipi.hpp')
f = h5py.File('pipi/pt_llll.1500.h5','r')
for i,ct in enumerate(f['pipi']['s0']):
  print(ct)

print('s1 in pipi.hpp')
f = h5py.File('pipi/pt_llll.1500.h5','r')
for i,ct in enumerate(f['pipi']['s1']):
  print(ct)


print('pipi correlator')
f = h5py.File('pipi/pt_llll.1500.h5','r')
for i,ct in enumerate(f['pipi']['corr']):
  print('d0={},  d1={},   c(t)=-d0+d1={}'.format(f['pipi']['d0'][i],f['pipi']['d1'][i],ct))
print()

print('dTst in pipi.hpp')
f = h5py.File('pipi/pt_llll.1500.h5','r')
for i,ct in enumerate(f['pipi']['cTst']):
  print(ct)


#print('s1 in pipi.hpp')
#f = h5py.File('pipi/pt_llll.1500.h5','r')
#for i,ct in enumerate(f['pipi']['cTst']):
#  print(ct)
#print('pipi correlator - pre-sinked')
#f = h5py.File('pipi/pt_smear_llll.1500.h5','r')
#for ct in f['pipi']['corr']:
#  print(ct)
#print()
