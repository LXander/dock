import os,sys
import pandas as pd
import numpy as np

def calculate_mean(filePath,k=10):
    '''
    columns in the file filename,inner_id,energy
    :param filePath:
    :param k:
    :return:
    '''

    print "convert mean ",filePath
    filename = os.path.basename(filePath)
    score_file = pd.read_csv(filePath)
    cutoff_position = score_file[score_file['inner_id'] <= k]

    gou = cutoff_position.groupby('filename')
    records = []
    for name, group in gou:
        records.append([name,np.mean(group['energy'])])

    reduced_score = pd.DataFrame(data = records,columns = ['filename','energy'])
    reduced_score['label'] = reduced_score['filename'].apply(lambda x:1 if x.find('actives')>=0 else 0)
    sorted_score = reduced_score.sort('energy')
    labels = sorted_score['label']


    points = []
    total = len(labels)
    ones = sum(labels)
    zeros = total - ones
    for i in range(total):
        points.append([1.0 * sum(labels[:i]) / ones, 1.0 * (i - sum(labels[:i])) / zeros])

    np_points = np.array(points)
    auc = np.trapz(np_points[:, 0], x=np_points[:, 1])

    with open(FLAGS.auc_log_file, 'a') as fout:
        fout.write(filename.split('.')[0] + '\t' + str(auc) + '\n')

    df = pd.DataFrame(data=np_points, columns=['tpr', 'fpr'])
    df.to_csv(os.path.join(FLAGS.destPath, filename.replace('.', 'mean_{}_curve.'.format(k))))


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
    labels = list(sorted_position['label'])

    points = []
    total = len(labels)
    ones = sum(labels)
    zeros = total - ones
    for i in range(total):
        points.append([1.0 * sum(labels[:i]) / ones, 1.0 * (i - sum(labels[:i])) / zeros])

    np_points = np.array(points)
    auc = np.trapz(np_points[:, 0], x=np_points[:, 1])

    with open(FLAGS.auc_log_file,'a') as fout:
        fout.write(filename.split('.')[0]+ '\t'+str(auc)+'\n')

    df = pd.DataFrame(data = np_points,columns=['tpr','fpr'])
    df.to_csv(os.path.join(FLAGS.destPath,filename.replace('.','_curve.')))



def run(index=None):
    fileNameList = os.listdir(FLAGS.sourcePath)
    print "index ",index
    if not index==None:
        print "index ",index
        calculate_mean(os.path.join(FLAGS.sourcePath,fileNameList[index]))
    else: 
        map(lambda filename:calculate_mean(os.path.join(FLAGS.sourcePath,filename)),fileNameList)

class FLAGS:
    sourcePath = '/n/scratch2/xl198/dude/frame/affinity'
    destPath = '/n/scratch2/xl198/dude/frame/mean_curve'
    auc_log_file = os.path.join(destPath,'mean_auc.txt')

if __name__ == '__main__':
    args = sys.argv
    print args 
    if len(args)>1:
        index =int(args[1])
        print "index ",index
    	run(index)
    else:    
        run()

