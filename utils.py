import warnings
import os, shutil
import numpy as np
import pandas as pd
from sklearn import cluster
from sklearn import manifold
from sklearn import decomposition
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs	
from scipy.cluster.hierarchy import dendrogram, linkage
warnings.filterwarnings('ignore', category=Warning)

def read_file(file, sep, header, out_dir, out_file_name, index_col='0'):
	index_col = int(index_col)
	if header == 'False':
		frame = pd.read_csv(file, sep=sep, header=None, index_col=index_col).iloc[:,index_col:]
	else:
		frame = pd.read_csv(file, sep=sep, index_col=index_col).iloc[:,index_col:]
	frame.to_hdf(out_dir + out_file_name + '.h5.tmp', 'slice_0')

def read_file_list(file_list, sep, header, out_dir, out_file_name, index_col='0'):
	for i in range(len(file_list)):
		read_file(file_list[i], sep, header, index_col, i, out_dir, out_file_name)

def slice_file(df_file, out_dir, out_file_name, slice_size='1000'):
	slice_size = int(slice_size	)
	df = pd.HDFStore(df_file)['slice_0']
	b = int(np.ceil(float(df.shape[0]) / float(slice_size)))
	digit = int(np.floor(np.log10(b)) + 1)
	for i in range(b):
		slice_name = str(i).zfill(digit)
		start = i*slice_size
		end = np.min([(i+1)*slice_size, df.shape[0]])
		slice_ = pd.DataFrame(data=df.iloc[start:end,:].values, index=df.index[start:end], columns=df.columns)
		slice_.to_hdf(out_dir + out_file_name + '.sliced.h5', 'slice_' + slice_name)
	pd.DataFrame(data=np.array(df.shape + (b,)), index=['row', 'col','slice']).to_hdf(out_dir + out_file_name + '.sliced.h5', 'attr')
	pd.DataFrame(data=df.columns, index=df.columns).to_hdf(out_dir + out_file_name + '.sliced.h5', 'cols')
	pd.DataFrame(data=df.index, index=df.index).to_hdf(out_dir + out_file_name + '.sliced.h5', 'rows')

def patch_file(df_file, out_dir, out_file_name):
	df = pd.HDFStore(df_file)['slice_0']
	df.to_hdf(out_dir + out_file_name + '.whole.h5', 'slice_0')
	pd.DataFrame(data=np.array(df.shape + (1,)), index=['row','col', 'slice']).to_hdf(out_dir + out_file_name + '.whole.h5', 'attr')
	pd.DataFrame(data=df.columns, index=df.columns).to_hdf(out_dir + out_file_name + '.whole.h5', 'cols')
	pd.DataFrame(data=df.index, index=df.index).to_hdf(out_dir + out_file_name + '.whole.h5', 'rows')

def calc_mi(arr1, arr2, bins, m):
	fq = np.histogram2d(arr1, arr2, bins=(bins, bins))[0]/float(m)
	sm = np.sum(fq*float(m), axis=1)
	tm = np.sum(fq*float(m), axis=0)
	sm = np.asmatrix(sm/float(sm.sum()))
	tm = np.asmatrix(tm/float(tm.sum()))
	sm_tm = np.matmul(np.transpose(sm), tm)
	div = np.divide(fq, sm_tm, where = sm_tm != 0, out = np.zeros_like(fq))
	ent = np.log(div, where = div != 0, out = np.zeros_like(div))
	agg = np.multiply(fq, ent, out = np.zeros_like(fq), where = fq != 0)
	return agg.sum()

def calc_mi_mat(mat1, mat2, bins, m, key, out_dir, out_file_name):
	df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
	for c in mat2.index:
		df.loc[mat1.index, c] = mat1.apply(calc_mi, axis=1, args=(mat2.loc[c,:], bins, m))
	df.to_hdf(out_dir + out_file_name + '.h5.tmp', str(key))

def merge_mi_mats(in_dir, in_common_name, n_slice, out_dir, out_file_name):
	mi_slices = [f for f in os.listdir(in_dir) if f.find(in_common_name + '_mi_') >= 0]
	n = ((n_slice + 1) * n_slice) / 2
	if(len(mi_slices) < n):
		print('[EROR] --> [MERG] Missing MI file(s).')
		exit()
	digit = int(np.floor(np.log10(n)) + 1)
	rows = []
	for i in range(n_slice):
		row_cols = []
		for j in range(i):
			idx = int(j * n_slice + i - (j * (j + 1)) / 2)
			key = 'mi_' + str(idx).zfill(digit)
			row_cols.append(pd.HDFStore(in_dir + in_common_name + '_' + key + '.h5.tmp')[key].T)
		for j in range(i, n_slice):
			idx = int(i * n_slice + j - (i * (i + 1)) / 2)
			key = 'mi_' + str(idx).zfill(digit)
			row_cols.append(pd.HDFStore(in_dir + in_common_name + '_' + key + '.h5.tmp')[key])
		rows.append(pd.concat(row_cols, axis=1))
	df = pd.concat(rows, axis=0)
	df.to_hdf(out_dir + out_file_name + '_mi.h5', 'mi')
	for f in mi_slices:
		os.remove(in_dir + f)

def norm_mi_mat(in_mat_file, out_dir, out_file_name):
	hdf = pd.HDFStore(in_mat_file)
	df = hdf['mi']
	diag = np.asmatrix(np.diag(df))
	if in_mat_file == out_dir + out_file_name + '_mi.h5':
		hdf.put('norm_mi', df / np.sqrt(np.matmul(diag.T, diag)))
	else:
		df.to_hdf(out_dir + out_file_name + '_mi.h5', 'mi')
		(df / np.sqrt(np.matmul(diag.T, diag))).to_hdf(out_dir + out_file_name + '_mi.h5', 'norm_mi')
	hdf.close()

def heatmap(fig_num, df, linkage_, matrix, out_dir, out_file_name, tag, dendro='True'):
	fig = plt.figure(fig_num, figsize=(16,16), dpi=300)
	grid = gs.GridSpec(1, 1, wspace=0.2, hspace=0.2)
	subgrid = gs.GridSpecFromSubplotSpec(2, 2, subplot_spec=grid[0], wspace=0.01, hspace=0.01, height_ratios=[1, 10], width_ratios=[10,1])
	if dendro == 'True':
		ax = plt.Subplot(fig, subgrid[0])
		dendrogram(linkage_, ax=ax, orientation='top', labels=df.index, leaf_font_size=2, color_threshold=0)
		ax.axis('off')
		fig.add_subplot(ax)
	ax1 = plt.Subplot(fig, subgrid[2])
	cax = ax1.matshow(matrix, cmap='YlOrRd', interpolation=None, aspect='auto')
	ax1.set_xticks([])
	ax1.set_yticks([])
	fig.add_subplot(ax1)
	ax3 = plt.Subplot(fig, subgrid[3])
	ax3.axis('off')
	cbar = fig.colorbar(cax, ax=ax3, shrink=1.5, aspect=20, fraction=0.5, pad=0)
	plt.savefig(out_dir + out_file_name + '_' + tag + '.pdf', bbox_inches='tight')

def heatmap2(fig_num, data, labels, out_dir, out_file_name, tag):
	if data is None:
		return
	fig = plt.figure(fig_num, figsize=(16,16), dpi=300)
	grid = gs.GridSpec(1, 1, wspace=0.2, hspace=0.2)
	subgrid = gs.GridSpecFromSubplotSpec(2, 2, subplot_spec=grid[0], wspace=0.001, hspace=0.1, height_ratios=[1, 20], width_ratios=[20, 1])
	sorted_labels = labels.sort_values()
	bar = np.vstack((np.array(sorted_labels) + 1, ) * 2) / np.max(np.unique(labels) + 1)
	ax1 = plt.Subplot(fig, subgrid[0])
	ax1.matshow(bar, cmap='jet', interpolation=None, aspect='auto', vmin=0, vmax=1)
	ax1.set_xticks([])
	ax1.set_yticks([])
	fig.add_subplot(ax1)
	ax2 = plt.Subplot(fig, subgrid[2])
	cax = ax2.matshow(data.loc[sorted_labels.index, sorted_labels.index], cmap='YlOrRd', interpolation=None, aspect='auto')
	ax2.set_xticks([])
	ax2.set_yticks([])
	fig.add_subplot(ax2)
	ax3 = plt.Subplot(fig, subgrid[3])
	ax3.axis('off')
	cbar = fig.colorbar(cax, ax=ax3, shrink=1.5, aspect=20, fraction=0.5, pad=0)
	plt.savefig(out_dir + out_file_name + '_' + tag + '.pdf', bbox_inches='tight')

def hclust(fig_num, in_mat_file, out_dir, out_file_name, plot='True'):
	hdf = pd.HDFStore(in_mat_file)
	df = 1 - hdf['norm_mi']
	n = df.shape[0]
	linkage_ = linkage(df, method='ward')
	leaves = dendrogram(linkage_, no_plot=True, orientation='top', labels=df.index, leaf_font_size=2, color_threshold=0)['leaves']
	leaves = np.asarray(leaves)
	leaves = df.index[leaves[leaves < n]]
	hc_matrix = df.loc[leaves, leaves]
	hc_matrix.to_hdf(out_dir + out_file_name + '_clust.h5', 'hclust')
	if plot == 'True':
		heatmap(fig_num, df, linkage_, hc_matrix, out_dir, out_file_name, 'hclust', 'True')
	hdf.close()

def scatter(fig_num, Y_tsne, out_dir, out_file_name, tag, facecolor='none', edgecolor='r', marker='o', marker_size=20):
	fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	plt.scatter(Y_tsne[:,0], Y_tsne[:,1], facecolor=facecolor, edgecolor=edgecolor, marker=marker, s=marker_size)
	plt.ylabel('MICA-2')
	plt.xlabel('MICA-1')
	plt.xticks([])
	plt.yticks([])
	plt.savefig(out_dir + out_file_name + '_' + tag + '.pdf',  bbox_inches='tight')

def scatter2(fig_num, data, out_dir, out_file_name, marker_size=20, marker='o'):
	if data is None:
		return
	fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	lab = np.unique(data.loc[:, 'label'])
	colors = plt.cm.jet(np.linspace(0, 1, len(lab)+1))
	for z in lab:
		X = data.loc[data.loc[:, 'label'] == z, ['X', 'Y']]
		plt.scatter(X.loc[:, 'X'], X.loc[:, 'Y'], facecolor=colors[z+1], s=marker_size, marker=marker, vmin=0, vmax=len(lab), label=str(z+1) + '(' + str(X.shape[0]) + ')', alpha=0.7)
		center = np.mean(X, axis=0)
		plt.scatter(center.loc['X'], center.loc['Y'], marker='o', c='white', alpha=0.7, s=100, edgecolor='k')
		plt.scatter(center.loc['X'], center.loc['Y'], marker='$%d$' % (z+1), c='black', alpha=0.7, s=80, edgecolor='k')
	plt.ylabel('MICA-2')
	plt.xlabel('MICA-1')
	plt.xticks([])
	plt.yticks([])
	plt.legend(loc='center left', bbox_to_anchor=(1,0.5), title='Clusters(' + str(data.shape[0]) + ')')
	plt.savefig(out_dir + out_file_name + '_clust_k' + str(len(lab)) + '.pdf', bbox_inches='tight')

def tsne(fig_num, data, max_dim, out_dir, out_file_name, tag, key, perplexity=30, plot='True'):
	tsne_ = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10, early_exaggeration=12.0)
	Y_tsne = tsne_.fit_transform(data.iloc[:,0:max_dim])
	pd.DataFrame(data=Y_tsne, index=data.index, columns=['X','Y']).to_hdf(out_dir + out_file_name + '_clust.h5', key)
	if plot == 'True':
		scatter(fig_num, Y_tsne, out_dir, out_file_name, tag)

def mds(fig_num, in_mat_file, max_dim, out_dir, out_file_name, plot='True'):
	hdf = pd.HDFStore(in_mat_file)
	df = 1 - hdf['norm_mi']
	n = df.shape[0]
	H = np.eye(n) - np.ones((n,n))/n
	B = -H.dot(df ** 2).dot(H)/2
	evals, evecs = np.linalg.eigh(B)
	idx = np.argsort(evals)[::-1]
	evals = evals[idx]
	evecs = evecs[:,idx]
	evals_pos = evals > 0
	L = np.diag(np.sqrt(evals[evals_pos]))
	V = evecs[:, evals_pos]
	LL = L[0:np.min([L.shape[0], 200]), 0:np.min([L.shape[0], 200])]
	VV = V[:, 0:LL.shape[0]]
	Y = pd.DataFrame(data=VV.dot(LL), index=df.index, columns=['mds_' + str(x) for x in np.arange(1, LL.shape[0] + 1)])
	Y.to_hdf(out_dir + out_file_name + '_clust.h5', 'mds')
	if plot == 'True':
		tsne(fig_num, Y, max_dim, out_dir, out_file_name, 'mds', 'mds-tsne')
	hdf.close()

def lpl(fig_num, in_mat_file,  max_dim, out_dir, out_file_name, plot='True'):
	hdf = pd.HDFStore(in_mat_file)
	df = 1 - hdf['norm_mi']
	n = np.min(df.shape[0], 200)
	laplacian = manifold.SpectralEmbedding(n_components=n, eigen_solver='lobpcg', random_state=10)
	Y = pd.DataFrame(data=laplacian.fit_transform(df), index=df.index, columns=['lpl_' + str(x) for x in np.arange(1, n)])
	Y.to_hdf(out_dir + out_file_name + '_clust.h5', 'lpl')
	if plot == 'True':
		tsne(fig_num, Y, max_dim, out_dir, out_file_name, 'lpl', 'lpl-tsne')
	hdf.close()

def pca(fig_num, in_mat_file, max_dim, out_dir, out_file_name, plot='True'):
	hdf = pd.HDFStore(in_mat_file)
	df = 1 - hdf['norm_mi']
	n = np.min(df.shape[0], 200)
	pca_ = decomposition.PCA(n_components=n, random_state=10)
	Y = pd.DataFrame(data=np.transpose(pca_.fit(df).components_), index=df.index, columns=['pca_' + str(x) for x in np.arange(1, n+1)])
	Y.to_hdf(out_dir + out_file_name + '_clust.h5', 'pca')
	if plot == 'True':
		tsne(fig_num, Y, max_dim, out_dir, out_file_name, 'pca', 'pca-tsne')
	hdf.close()

def kmeans(in_mat_file, transformation, n_clusters, dim, out_dir, out_file_name):
	hdf = pd.HDFStore(in_mat_file)
	df = hdf[transformation]
	km = cluster.KMeans(n_clusters=n_clusters, max_iter=1000, n_init=1000)
	km_res = pd.DataFrame(data=np.transpose(km.fit_predict(df.iloc[:,0:dim])), index=df.index, columns=['label'])
	km_res.to_hdf(out_dir + out_file_name, 'kmeans')
	hdf.close()

def aggregate(tmp_dir, n_clusters, common_name):
	clusts = [f for f in os.listdir(tmp_dir) if f.find(common_name) >= 0]
	n_iter = len(clusts)
	if n_iter == 0:
		return None
	mem = None
	for i in range(n_iter):
		hdf = pd.HDFStore(tmp_dir + clusts[i]) 
		df = hdf['kmeans'] + 1
		dff = pd.DataFrame(data=np.matmul(df, df.T), index=df.index, columns=df.index)
		dff_div = pd.DataFrame(data=np.array((np.diag(dff),)*dff.shape[0]).T, index=dff.index, columns=dff.columns)
		mem_mat = pd.DataFrame(data=dff/dff_div==1, index=dff.index, columns=dff.columns, dtype=np.float32)
		mem = mem_mat if i == 0 else mem + mem_mat.loc[mem.index, mem.columns]
		mem = mem / n_iter if i == n_iter - 1 else mem
		hdf.close()
	clust = cluster.AgglomerativeClustering(linkage='ward', n_clusters=n_clusters, affinity='euclidean')
	clust.fit(mem)
	cclust = pd.DataFrame(data=clust.labels_, index=mem.index, columns=['label'])
	index = cclust.groupby(['label']).size().sort_values(ascending=False).index
	label_map = {index[i]: i+1000 for i in range(len(index))}
	cclust = cclust.replace(to_replace={'label': label_map})
	index = cclust.groupby(['label']).size().sort_values(ascending=False).index
	map_back = {index[i]: i for i in range(len(index))}
	cclust = cclust.replace(to_replace={'label': map_back})
	return [cclust, mem]

def cc(tmp_dir, n_clusters, project_name, common_name, transformation, out_dir, out_file_name, max_dim=0, re_trans='False'):
	agg = aggregate(tmp_dir, n_clusters, common_name)
	if agg != None:
		cclust = agg[0]
		mem = agg[1]
		out_file = out_dir + out_file_name + '_clust.h5'
		if os.path.exists(out_file):
			transformation = 'pca' if transformation == 'lpca' else transformation
			if re_trans == 'True' and max_dim > 0:
				hdf = pd.HDFStore(out_file)
				Y = hdf[transformation]
				hdf.close()
				perplexity = np.min([1200, np.max(cclust.groupby(['label']).size())])
				tsne(None, Y, max_dim, out_dir, out_file_name, None, transformation + '-tsne', perplexity, 'False')
				hdf = pd.HDFStore(out_file)
				Y_tsne = hdf[transformation + '-tsne']
				hdf.close()
			else:
				hdf = pd.HDFStore(out_file)
				Y_tsne = hdf[transformation + '-tsne']
				hdf.close()
			cclust = pd.concat([cclust, Y_tsne.loc[cclust.index, :]], axis=1)
			data_file = tmp_dir + project_name + '.whole.h5'
			if os.path.exists(data_file):
				hdf = pd.HDFStore(data_file)
				data = hdf['slice_0']
				cclust = pd.concat([cclust, data.loc[cclust.index, :]], axis=1)
				hdf.close()
			mi_file = tmp_dir + project_name + '_mi.h5'
			if os.path.exists(mi_file):
				hdf = pd.HDFStore(mi_file)
				mi = 1 - hdf['norm_mi']
				mi.to_hdf(out_file, 'norm_mi')
				hdf.close()
		cclust.to_hdf(out_file, 'cclust')
		mem.to_hdf(out_file, 'membership')
		shutil.rmtree(tmp_dir)

def plot(fig_num, in_file, out_dir, out_file_name, hclust='False'):
	hdf = pd.HDFStore(in_file)
	mi = hdf['norm_mi'] if '/norm_mi' in hdf.keys() else None
	mem = hdf['membership'] if '/membership' in hdf.keys() else None
	cclust = hdf['cclust'] if '/cclust' in hdf.keys() else None
	scatter2(fig_num, cclust, out_dir, out_file_name)
	if cclust is not None and hclust == 'True':
		heatmap2(fig_num + 1, mi, cclust.loc[:, 'label'], out_dir, out_file_name, 'mi')
		heatmap2(fig_num + 2, mem, cclust.loc[:, 'label'], out_dir, out_file_name, 'clusts')
	hdf.close()


