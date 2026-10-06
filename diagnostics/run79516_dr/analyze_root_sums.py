import uproot,numpy as np,json,os
for tag in ['true','false']:
 path='/tmp/cpm_'+tag+'.root'
 if not os.path.exists(path):continue
 f=uproot.open(path);a=f['cpm_voxel_correction_sums'].arrays(library='np');n=a['entries'];dr=a['sum_delta_r']/n;dca=a['sum_dca']/n;r=np.hypot(a['sum_voxel_x']/n,a['sum_voxel_y']/n)
 print('\nFILE',tag,'voxels',len(n),'pairs',int(n.sum()))
 print('summary',{k:v.tolist() for k,v in f['cpm_b3_summary'].arrays(library='np').items() if k in ['max_pair_dca','crossing_solver','magnetic_field_z','accepted_pairs','candidate_pairs']})
 print('dr quantiles',np.quantile(dr,[0,.01,.1,.5,.9,.99,1]),'weighted mean',np.sum(a['sum_delta_r'])/n.sum())
 print('dca quantiles',np.quantile(dca,[0,.5,.99,1]))
 print('voxel fraction abs dr>10',np.mean(abs(dr)>10),'pairs fraction in those',n[abs(dr)>10].sum()/n.sum())
 h=f['hDistortionR_rec'].values();print('hist vs sums max diff',np.max(abs(h[a['iphi'],a['ir'],a['iz']]-dr)))
 print('plain weighted match',np.max(abs(a['sum_weighted_delta_r']/a['sum_pair_weight']-dr)))
 for j in np.argsort(abs(dr))[-5:][::-1]:print('extreme', {k:float(v[j]) for k,v in {'r':r,'dr':dr,'mean_midpoint_r':r-dr,'n':n,'dca':dca,'rms_dr':np.sqrt(np.maximum(0,a['sum_delta_r2']/n-dr**2)),'iphi':a['iphi'],'ir':a['ir'],'iz':a['iz']}.items()})
 print('radial rows')
 for ir in np.unique(a['ir']):
  m=a['ir']==ir;print(int(ir),round(np.average(r[m],weights=n[m]),2),round(np.sum(a['sum_delta_r'][m])/n[m].sum(),3),int(n[m].sum()))
