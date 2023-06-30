def ChangePointDetection():
    import os
    import scipy.io as sio
    import scipy.stats as stats
    import numpy as np
    import ruptures as rpt
    import matplotlib.pylab as plt
    from IPython import display
    import time
    
    os.chdir(r"C:\Users\dmagyar\Desktop\M2_shock_waveform_ChangePointDtct\kilosort\MD129_kilosort\kilosort3preprocess")
    mat = sio.loadmat('temp_wh.cell_metrics.cellinfo.mat')
    cell_metrics = np.array(mat['cell_metrics'])
    manipulations = cell_metrics[['manipulations']]
    cellID = cell_metrics[['cellID']]
    for ii in range(cellID[0][0][0].size):
        curr_manip = manipulations[0][0][0][0][0][0][0][ii]
        curr_cellID = cellID[0][0][0][0][ii]
        z_scores = stats.zscore(curr_manip)
        offset = int((z_scores.size-1)/4) # hardcoded, psth baseline used in CellExplorer
        model="rank"
        std_noise = np.std(z_scores[:offset])
        std_noise_full =np.std(z_scores)
        pen = np.log(z_scores.size) * 1 * std_noise**2
        algo = rpt.Pelt(model=model,jump=1).fit(np.array(z_scores))
        result = algo.predict(pen=pen)
        f_result = [x for x in result if (150 > x > 45 )][:4]
        pointchck = "unknown"
        while pointchck != "++":
            currPlot = rpt.display(np.array(z_scores), f_result, figsize=(10, 6))
            currPlot[0].canvas.draw()
            plt.ion()
            plt.show(block=False)
            
            print(f"Currently keeping changepoints:{f_result}")
            print(f"Remove: (first index is 0)")
            pointchck = input()
            if pointchck != "++":
                pointchck = int(pointchck)
                f_result.pop(pointchck)
                plt.close()
            

ChangePointDetection()
