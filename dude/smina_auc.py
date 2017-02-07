import os,sys
import pandas as pd
import numpy as np

def calculate(filePath):
    '''
    Given csv store smina score, calculate the auc based on the best position's smina score
    :param filePath:
    :return:
    '''

    print "convert ",filePath
    filename = os.path.basename(filePath)
    score_file = pd.read_csv(filePath)
    first_position = score_file[score_file['inner_id']==1]
    sorted_position = first_position.sort('energy')
    sorted_position['label'] = sorted_position['filename'].apply(lambda x:1 if x.find('actives')>=0 else 0)
    labels = list(sorted_res['label'])

    points = []
    total = len(labels)
    ones = sum(labels)
    zeros = total - ones
    for i in range(total):
        points.append([1.0 * sum(labels[:i]) / ones, 1.0 * (i - sum(labels[:i])) / zeros])

    np_points = np.array(points)
    auc = np.trapz(np_points[:, 0], x=np_points[:, 1])

    with open(FLAGS.auc_log_file,'a') as fout:
        fout.write(filename.split('.')+ '\t'+str(auc)+'\n')

    df = pd.DataFrame(data = np.points,columns=['tpr','fpr'])
    df.to_csv(os.path.join(FLAGS.destPath,filename.replace('.','_curve.')))



def run():
    fileNameList = os.listdir(FLAGS.sourcePath)
    map(lambda filename:calculate(os.path.join(FLAGS.sourcePath,filename)),fileNameList)

class FLAGS:
    sourcePath = '/n/scratch2/xl198/dude/frame/affinity'
    auc_log_file = '/n/scratch2/xl198/dude/frame/curve/auc.txt'
    destPath = '/n/scratch2/xl198/dude/frame/curve'

if __name__ == '__main__':
    run

