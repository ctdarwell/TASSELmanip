#input file generated in R using: mkHM.R
    
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cv2, itertools
from scipy import ndimage

#PARAMS - can be set from vcf2ldheatmap_FromALL.py
write = [True, False][0]
file = 'chrom_11.csv'
workfldr = '.'
prefix = 'TEST'
dpi = 300
region = [17000000, 29000000] # = [start, end] OR None; chr region shown on img
genes = ['LOC_Os11g31620', 'LOC_Os11g32210', 'Xa21', 'LOC_Os11g36390', 'Xa10/Xa39', 
         'LOC_Os11g37740', 'LOC_Os11g37860', 'LOC_Os11g37870', 'LOC_Os11g37960', 
         'LOC_Os11g38870', 'LOC_Os11g40840', 'LOC_Os11g45620', 'Xa40/Xa3/Xa26'] #match to 'markers_to_highlight' positions
markers_to_highlight = [18486667, 19025182, 20805000, 21446864, 22260000, 
                      22311354, 22432467, 22444386, 22507045, 
                      23136993, 24438215, 27605223, 28301000] #add more according to 'genes'
colours = ['dkred', 'dkred', 'black', 'dkred', 'black', 
           'dkred', 'dkred', 'dkred',  'dkred',
            'dkred', 'dkred', 'dkred','black']
rois = [[17438481, 19648392], [20894905, 24947194]] # set regions of interest [list of lists] - draws red triangles on heatmap

colDik = {'white': (255, 255, 255), 'blue': (200, 0, 0), 'red': (0, 0, 255), 'black': (1, 0, 0), 'green': (0, 200, 0), 'dkred': (0, 0, 150)}
roicol = colDik.get('dkred') #can change region of interest colour
pointercols = []
[pointercols.append(colDik.get(c)) for c in colours]
pointerdik = dict(zip(markers_to_highlight, pointercols))

def main():
    '''
    Main pgm control function
    '''
    if (len(genes) != len(markers_to_highlight)) or (len(genes) != len(pointercols)) or (len(markers_to_highlight) != len(pointercols)):
        print('RUN STOPPED: variables "genes", "markers_to_highlight", and "pointercols" are not all of equal length')
        return
    pointers = markers_to_highlight
    dat = pd.read_csv(f"{workfldr}/{file}", header=0)
    print(f"Read: {file}; size = {dat.shape}")
    im, cols = mkHMP(dat)
    starty, startx, endx = getXYcoords(im)
    snpmap = [[startx, endx], region]
    maplen = snpmap[1][1] - snpmap[1][0]
    loci = dat.columns[cols[:-1]]
    rois, pointers = checkMarkers(loci, pointers)
    im, pointerlen, pointerWid, pointerLocs = mkPointers(im, pointerdik, snpmap, maplen, loci, startx, endx, starty, pointers)
    im = mkROIs(im, starty, loci, startx, endx, rois)
    pointerLocs = unSquishTxt(im, pointerLocs, snpmap) #MIGHT GIVE PROBS AT SOME POINT!!
    im = writeText(im, pointerLocs, pointerlen, pointerWid)
    if write: cv2.imwrite(f"{workfldr}/{prefix}_{int(region[0]/1000000)}-{int(region[1]/1000000)}Mb_heatmap.png", im)
    print("KAT-AM!")

def mkHMP(dat):
    '''make the heat map'''
    print('\nHeatmap to be annotated [MAY BE BLANK!!]')
    cols = []
    if region:
        for col in dat.columns[1:]: #id df columns in chosen region (also useful later on)
            if float(col) > region[0]: cols.append(dat.columns.tolist().index(col))
            if float(col) > region[1]: break
    else: cols = [0, len(dat.columns) - 1]
    arr = np.array(dat.iloc[cols[0]:cols[-1], cols[0] + 2:cols[-1] + 2]) #extract stipulated region
    mask = np.zeros_like(arr)
    mask[np.tril_indices_from(mask)] = True #something to do with stating upper tri as lower
    
    #create a heatmap - cbar=True gives colorbar (fao extracting colour bar - heatmap is too small)
    fig = plt.figure(1)
    with sns.axes_style("white"):
        sns.heatmap(arr, xticklabels=False, cbar=True, 
                    cbar_kws={"orientation": "horizontal", "shrink": 0.18},
                    yticklabels=False, mask=mask, square=True,  cmap="YlOrRd")

    fig.savefig("tmp.png", dpi = dpi) #save so can open with cv2
    plt.show()
    plt.close(fig)

    #extract color bar
    im = cv2.imread('tmp.png')
    bar, ys, xs = getBarCoords(im)
    
    #create 2nd heatmap - cbar=False gives no colorbar (the heatmap is now bigger)
    fig = plt.figure(1)
    with sns.axes_style("white"):
        sns.heatmap(arr, xticklabels=False, cbar=False, 
                    yticklabels=False, mask=mask, square=True,  cmap="YlOrRd")

    fig.savefig("tmp.png", dpi = dpi)
    plt.close(fig)

    #rotate heatmap and add colour bar
    im = ndimage.rotate(cv2.imread("tmp.png"), 225, mode='nearest') #rotate by 225 deg
    y, x = int(im.shape[0] * .75), int(im.shape[0] * .565) #calc where to put bar
    im[y:y+bar.shape[0], x:x+bar.shape[1], :] = bar #overlay bar

    return im, cols

def getBarCoords(im2): #this could feasibly be mgd w getXYcoords as a single func
    '''similar to "getXYcoords" but finds outer edges of the colour bar
    then returns its coords'''
    l1 = []
    for starty2 in range(im2.shape[0] - 1, 0, -1): #pixel by pixel up image's centre til non-white pixel found 
        if im2[starty2, im2.shape[1]//2, 0] != 255 or im2[starty2, im2.shape[1]//2, 1] != 255 or im2[starty2, im2.shape[1]//2, 2] != 255:
            l1.append(starty2)

    for l in l1: #evaluate top & bottom border of bar - list of non-white pixels need splitting into bar vs. heatmap
        if l - l1[l1.index(l) + 1] > im2.shape[0] / 8: break #finds where the bar ends (going upwards) - i.e. the gap btwn the bar and the actual heat map (there are gaps of white within the bar so this is a bit of fiddling to find the real gap)
    ys = [l1[0] + int(im2.shape[0] / 80), l - int(im2.shape[0] / 80)] #top & bottom of bar

    for startx2 in range(im2.shape[1]): #find LH edge of bar using upper edge y-coord of bar
        if im2[l, startx2, 0] != 255 or im2[l, startx2, 1] != 255 or im2[l, startx2, 2] != 255:
            break

    xs = [startx2 - int(im2.shape[0] / 10), im2.shape[1] - startx2 + int(im2.shape[0] / 10)]
    bar = im2[ys[1]: ys[0], xs[0]:xs[1]] #coords of bar with a border
    
    return bar, ys, xs

def getXYcoords(im):
    '''identify start/end posns (i.e. margines) of heat map against blank background
    work down y-axis from midpoint of x (we know heat map is in centre of img
    from id'd y, work inwards to find LHS x value
    endx given by symmetry of shape'''
    for starty in range(im.shape[0]): #pixel by pixel down image's centre til non-white pixel found
        if im[starty, im.shape[1]//2, 0] != 255 or im[starty, im.shape[1]//2, 1] != 255 or im[starty, im.shape[1]//2, 2] != 255:
            break

    for startx in range(im.shape[1]): #find LH edge of heatmap using upper edge y-coord of bar
        if im[starty, startx, 0] != 255 or im[starty, startx, 1] != 255 or im[starty, startx, 2] != 255:
            break

    endx = im.shape[1] - startx #RH edge of heatmap

    return starty, startx, endx

def checkMarkers(loci, pointers):
    '''check stipulated markers are actually in the data -
    if not, for ROIs, take markers giving widest region within the stipulated boundaries
    for pointers, take nearest marker - but keep record of original input as this chsomal posn will be mapped'''

    newrois, newpointers = [], []
    mks = loci.astype(np.int)
    for roi in rois:
        if roi[0] not in mks: print(f"No recorded SNP at {roi[0]}!\nAltering region of interest values accordingly.")
        if roi[1] not in mks: print(f"No recorded SNP at {roi[1]}!\nAltering region of interest values accordingly.")
        newrois.append([mks[np.where(mks >= roi[0])[0][0]], mks[np.where(mks > roi[1])[0][0] - 1]]) #get next inward SNPs from stipulated start-end of ROI
    print('\n')
    for pointer in pointers: 
        if pointer not in mks: #below will store two values of pointer as list within list of pointers (i.e. inputted value shows chromosomal posn while calculated nearest value used to mark on LD map) 
            print(f"No recored SNP at {pointer}!\nAltering gene position to nearest marker")
            v1, v2 = abs(pointer - mks[np.where(mks > pointer)[0][0] - 1]), abs(pointer - mks[np.where(mks > pointer)[0][0]])
            if v1 < v2: newpointers.append([pointer, mks[np.where(mks > pointer)[0][0] - 1]])
            else: newpointers.append([pointer, mks[np.where(mks > pointer)[0][0]]])
        else: newpointers.append(pointer)

    return newrois, newpointers

def mkPointers(im, pointerdik, snpmap, maplen, loci, startx, endx, starty, pointers):
    '''calculate start and end point of pointers and draw.
    Bottom of pointer is at a SNP featured in the data (maybe recalculated in 'checkMarkers' if non-present marker given).
    Top of pointer is original given position, irrespective of whether the position has a SNP in the data or not'''
    pointerLocs = []
    pointerlen, pointerGap, pointerWid = int(im.shape[0] / 20), int(im.shape[0] / 280), int(im.shape[0] / 350)
    for p in pointers:
        try: p0, p1 = p[0], p[1] #if pointer has recalculated value alongside original that has no representative marker
        except: p0, p1 = p, p
        pointprop  = (p0 - snpmap[1][0]) / maplen 
        pointerTOP = int(startx + (pointprop * (endx - startx))) #start point is proportional to end/start of map
        pointerBOT = startx + int(np.where(loci == str(p1))[0][0] / loci.shape[0] * (endx - startx)) #bottom is proportional according to no. of markers
        pointerLocs.append([pointerTOP, starty - pointerlen - pointerGap]) #need a record as referenc for later drawing
        incr = (pointerBOT - pointerTOP) / pointerlen #calc x increment per y pixel
        for ypos in range(starty - pointerlen - pointerGap, starty - pointerGap): #draw the lines
            im[ypos:ypos + 1, int(pointerTOP + incr): int(pointerTOP + incr + pointerWid)] = pointerdik.get(p0)
            pointerTOP += incr

    return im, pointerlen, pointerWid, pointerLocs

def mkROIs(im, starty, loci, startx, endx, rois):
    '''draw triangles marking regions of interest
    maybe recalculated in 'checkMarkers' if locations given don't have a SNP in the data'''
    for roi in rois: #again calc like pointerBOT
        roistart = startx + int(np.where(loci == str(roi[0]))[0][0] / loci.shape[0] * (endx - startx))
        roiend = startx + int(np.where(loci == str(roi[1]))[0][0] / loci.shape[0] * (endx - startx))
        lim1, lim2 = int(im.shape[0] / 212), int(im.shape[0] / 170)

        #draw a hollow triangle
        y2 = starty + 5
        l, r, d = roistart, roiend, 1
        while d < lim1: #width of horizontal line of tri
            im[y2 + d -1: y2 + d, l:r] = roicol
            l += 1
            r -= 1
            d += 1
        while r - lim2 > l: #draw angled lines of tri with gaps between (ie don't overwrite heatmap)
            im[y2 + d -1: y2 + d, l:l + lim2] = roicol
            im[y2 + d -1: y2 + d, r - lim2:r] = roicol
            l += 1
            r -= 1
            d += 1
        while r > l: #make tip w/o gap btwn lines
            im[y2 + d -1: y2 + d, l:r] = roicol
            l += 1
            r -= 1
            d += 1

    return im

def unSquishTxt(im, pointerLocs, snpmap):
    '''check there is enough space between the labels
    re-position if too squished'''
    gap = int(im.shape[0] / 50)
    dists, probs, cnt = [], [[]], 0
    [dists.append(x[0]) for x in pointerLocs]
    #dists.insert(0, snpmap[0][0]), dists.append(snpmap[0][1]) #OLD LINE
    dists.insert(0, 100), dists.append(im.shape[1] - 100) #ADDED EXPT LINE - imaginery locations near the limits of the im so calcs below have values to check if they're near something (i.e. the list indexes don't get out of bounds) - a bit fudgey
    for dist in dists[:-1]:
        #make list of problematic posns
        if dists[dists.index(dist) + 1] - dist < gap: probs[cnt].extend([dist, dists[dists.index(dist) + 1]])
        else:
            probs.append([])
            cnt += 1

    #house keeping
    while [] in probs: probs.remove([])
    if probs == []: return pointerLocs
    
    swch, cntr = [1], 0
    while swch != []:
        cntr += 1
        if cntr > 100:
            print('RUN STOPPED! There is a problem seperating labels on the figure.\nTry expanding the "region" variable to allow more space')
            return
        swch = []
        #merge any clusters of problems together if they contain a common problem posn
        for it in itertools.combinations(probs, 2):
            if list(set(it[0]) & set(it[1])) != []:
                probs.remove(probs[probs.index(it[1])])
                probs.remove(probs[probs.index(it[0])])
                probs.insert(0, it[0] + it[1])
                
        for prob in probs:
            #calculate new gaps btwn txt (i.e. x coords)
            setprob = sorted(list(set(prob)))
            mu = int(np.mean([setprob[0], setprob[-1]]))
            init = mu - ((len(setprob) - 1) * int(gap / 2))
            for i in setprob[:-1]: init += gap
            #check new posns aren't problematic - otherwise back to top of loop and redo w new posns
            if mu - ((len(setprob) - 1) * int(gap / 2)) - dists[dists.index(setprob[0]) - 1] < gap:
                probs[probs.index(prob)].append(dists[dists.index(setprob[0]) - 1])
                swch.append(1) #any issues mean swch != []: thus, goes back to top of while loop
            if dists[dists.index(setprob[-1]) + 1] - init < gap:
                probs[probs.index(prob)].append(dists[dists.index(setprob[-1]) + 1])
                swch.append(1)

    #assign new gaps
    for prob in probs:
        setprob = sorted(list(set(prob))) #recalc these
        mu = int(np.mean([setprob[0], setprob[-1]]))
        init = mu - ((len(setprob) - 1) * int(gap / 2))
        for i in setprob: #change vals in pointerLocs
            pointerLocs[pointerLocs.index([i, pointerLocs[0][1]])][0] = init
            init += gap

    return pointerLocs    

def writeText(im, pointerLocs, pointerlen, pointerWid):
    '''write text on dummy array same size as im
    then rotate before applying to im with appropriate colors'''
    xaposns = dict(zip(genes, pointerLocs))
    font = cv2.FONT_HERSHEY_SIMPLEX
    fontScale = int(im.shape[0] / 2000)
    lineType = int(im.shape[0] / 500)

    for gene in genes:
        txt = np.zeros(im.shape, np.uint8) # Create a black image - NB pointerlen determiness dist from pointer!!
        bottomLeftCornerOfText = (xaposns.get(gene)[1] + int(pointerlen * 2),
                                  xaposns.get(gene)[0] + pointerWid) #txtLoc(xaposns.get(gene))
        fontColor = pointercols[genes.index(gene)]
        cv2.putText(txt, gene, 
                    bottomLeftCornerOfText, 
                    font, 
                    fontScale,
                    fontColor,
                    lineType)

        txt = ndimage.rotate(txt, 90, mode='nearest')

        txtposns = np.where(txt != 0) #get txt posns in txt (i.e. that aren't white) and overwrite only these on im 

        im[txtposns[0], txtposns[1], :] = fontColor

    return im

if __name__ == '__main__': main()
