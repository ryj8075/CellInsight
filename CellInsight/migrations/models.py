import boto3
from botocore.client import Config

import io
import os

import gzip
import h5py
import anndata
import pandas as pd
import scanpy as sc

from rpy2 import robjects
from rpy2.robjects.packages import importr

bucket_name='cellinsight-bucket'

# boto3.resource provides low-level interface 

s3_resource = boto3.resource(
    's3',
    aws_access_key_id='ncp_iam_BPAMKR1DK9BTckMdUIOn',  
    aws_secret_access_key='ncp_iam_BPKMKRZSqttkGfV9d8hY6sZ5Hg057NU9E3',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)

# boto3.client provide high-levee interface(object-oriented interface)

s3_client = boto3.client(
    's3',
    aws_access_key_id='ncp_iam_BPAMKR1DK9BTckMdUIOn',  
    aws_secret_access_key='ncp_iam_BPKMKRZSqttkGfV9d8hY6sZ5Hg057NU9E3',  
    endpoint_url='https://kr.object.ncloudstorage.com',
    region_name = 'kr-standard',
    config=Config(signature_version='s3v4')
)

def list_s3_objects(bucket_name, prefix, delimiter=False):
    try:
        params = {
            'Bucket': bucket_name,
            'Prefix': prefix,
        }
        
        if delimiter:
            params['Delimiter'] = '/'

        response = s3_client.list_objects_v2(**params)
        return response
    except:
        return None

def get_s3_objects(bucket_name, file_key):
    try:
        response = s3_client.get_object(Bucket=bucket_name, Key=file_key)
        unprocessed_data = response['Body'].read()
        return unprocessed_data
    except:
        return None
    
def download_s3_objects(bucket_name, file_key, file_path):
    try:
        s3_client.download_file(bucket_name, file_key, file_path)
        return 
    except:
        return False

def read_RData(file_data):
    base = importr('base')
    
    with robjects.conversion.localconverter(robjects.default_converter + robjects.numpy2ri.converter):
        r_object = base.readRDS(io.BytesIO(file_data))
        return robjects.conversion.rpy2py(r_object)

def read_h5ad(file_data):
    with io.BytesIO(file_data) as f:
        data = anndata.read_h5ad(f)
        return data
    
def read_scp_cluster(file_data):
    try:
        # 만약 file_data가 bytes 객체라면, 이를 파일 객체로 변환
        if isinstance(file_data, bytes):
            file_data = io.BytesIO(file_data)
    
        # 파일의 첫 1024 바이트를 읽어서 파일 형식을 판별합니다.
        file_snippet = file_data.read(1024)  # 첫 1024 바이트 읽기 (바이트 타입)
        file_data.seek(0)  # 파일 포인터를 처음으로 되돌립니다.
    
        # 파일 형식 판별 및 처리
        if file_snippet.startswith(b'\x1f\x8b'):  # gzip 파일
            # gzip 파일 처리
            try:
                with gzip.open(file_data, 'rt', encoding='utf-8') as gz_file:
                    df = pd.read_csv(gz_file)
            except Exception as e:
                print(f"Error reading gzip file: {e}")
                return None
        
        elif file_snippet[:4] == b'\x89HDF':  # HDF5 파일 (일부 형식 예: H5 or H5AD)
            print("Cluster files cannot be HDF5 format.")
            return None  # Return None for unsupported HDF5 format
        
        else:
            # 텍스트 파일 처리
            file_snippet = file_snippet.decode('utf-8')  # 바이트를 문자열로 변환
        
            # 파일 확장자에 따라 적절한 구분자 설정
            if ',' in file_snippet:
                sep = ',' 
            else:
                sep = '\t'  # 기본적으로는 탭으로 구분된 파일로 가정

            try:
                # 파일을 데이터프레임으로 읽어옵니다.
                df = pd.read_csv(file_data, sep=sep)
            except Exception as e:
                print(f"Error reading CSV/TSV file: {e}")
                return None  # Return None if an error occurs

        return df
    
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def read_scp_exp(file_data):
    
    try:
        # 만약 file_data가 bytes 객체라면, 이를 파일 객체로 변환
        if isinstance(file_data, bytes):
            file_data = io.BytesIO(file_data)

        # 파일의 첫 4바이트를 읽어 파일 형식을 판별합니다.
        file_snippet = file_data.read(4)
        file_data.seek(0)  # 파일 포인터를 처음으로 되돌립니다.

        # h5 또는 h5ad 파일인 경우 처리 (HDF5 파일은 0x89HDF로 시작)
        if file_snippet == b'\x89HDF':
            try:
                with h5py.File(file_data, 'r') as f:
                    # H5 파일이 h5ad 또는 10X Genomics 형식인지 확인합니다.
                    if 'matrix' in f.keys() or 'X' in f.keys():
                        return sc.read_10x_h5(file_data)
                    else:
                        return sc.read_h5ad(file_data)
            except Exception as e:
                print(f"Error reading h5/h5ad file: {e}")
                return None 

        # 텍스트 파일을 처리할 경우 (csv/tsv)
        # 파일의 첫 1024 바이트를 다시 읽고 문자열로 변환하여 구분자를 판별합니다.
        file_snippet = file_data.read(1024).decode('utf-8')
        file_data.seek(0)  # 파일 포인터를 처음으로 되돌립니다.

        # 파일 확장자에 따라 적절한 구분자 설정
        if '\t' in file_snippet:  # tsv
            sep = '\t'
        elif ',' in file_snippet:  # csv
            sep = ',' 
        else:
            sep = '\t'  # 기본적으로는 탭으로 구분된 파일로 가정

        # 텍스트 파일을 읽기 위해 파일 데이터를 문자열로 디코딩하여 StringIO로 변환
        file_data_str = io.StringIO(file_data.read().decode('utf-8'))
        file_data.seek(0)  # 파일 포인터를 처음으로 되돌립니다.

        try:
            # 파일을 anndata로 읽어옵니다.
            adata = sc.read_csv(file_data_str, delimiter=sep)
            return adata
        except Exception as e:
            print(f"Error reading CSV/TSV file: {e}")
            return None
    
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def read_scp_exp_gz(file_data):
    with gzip.GzipFile(fileobj=io.BytesIO(file_data)) as gz:
        # gzip으로 압축된 파일을 읽고, 읽은 데이터를 그대로 read_scp_exp 함수로 전달
        decompressed_data = gz.read()
        return read_scp_exp(decompressed_data)
    
def read_10x_mtx(folder_path):
    # 폴더 안의 모든 파일을 탐색
    for filename in os.listdir(folder_path):
        old_file_path = os.path.join(folder_path, filename)
        
        # 파일이 아닌 디렉토리일 경우 건너뜀
        if not os.path.isfile(old_file_path):
            continue
        
        # 파일 확장자와 이름 검사를 위한 초기화
        new_filename = None
        
        # 확장자가 '.mtx' 또는 '.mtx.gz'인 경우
        if filename.endswith('.mtx'):
            new_filename = 'matrix.mtx'
        elif filename.endswith('.mtx.gz'):
            new_filename = 'matrix.mtx.gz'
        
        # 확장자가 '.tsv' 또는 '.tsv.gz'인 경우
        elif filename.endswith('.tsv'):
            if 'genes' in filename:
                new_filename = 'genes.tsv'
            else:
                new_filename = 'barcodes.tsv'
        elif filename.endswith('.tsv.gz'):
            if 'genes' in filename:
                new_filename = 'genes.tsv.gz'
            else:
                new_filename = 'barcodes.tsv.gz'
        
        # 파일 이름이 변경되어야 할 경우, 이름 변경 수행
        if new_filename:
            new_file_path = os.path.join(folder_path, new_filename)
            
            # 동일한 이름의 파일이 이미 존재할 경우, 파일 이름 변경을 건너뜀
            if not os.path.exists(new_file_path):
                os.rename(old_file_path, new_file_path)
                print(f"Renamed: {filename} -> {new_filename}")
            else:
                print(f"Skipped: {new_filename} already exists")
        else:
            print(f"No change needed for: {filename}")
    
    return sc.read_10x_mtx(folder_path, var_names='gene_symbols')
    
def read_PanglaoDB_from_s3(bucket_name, prefix):
    response = list_s3_objects(bucket_name, prefix)
    
    if response == None:
        return None
    
    file_list = [] 
    
    for file in response['Contents']:
        file_key = file['Key']
        unprocessed_data = get_s3_objects(bucket_name, file_key)
        
        if file_key.endswith('.RData'):
            processed_data = read_RData(unprocessed_data)

            parts = file_key.split('/')
            information = {'study' : parts[1], 'file_name' : parts[2], 'data': processed_data}
                
            file_list.append(information)
    
    return file_list

def process_exp(bucket_name, folder_prefix):
    response = list_s3_objects(bucket_name, folder_prefix)
        
    if response is None or 'Contents' not in response:
        # S3 객체가 없거나, 응답에 Contents가 없으면 빈 리스트 반환
        return []
        
    file_list = []
        
    for file in response['Contents']:
        file_key = file['Key']
        unprocessed_data = get_s3_objects(bucket_name, file_key)
        success = False
            
        # Check file type and process accordingly
        if file_key.endswith('.gz'):
            processed_data = read_scp_exp_gz(unprocessed_data)
            success = True
        elif any(file_key.endswith(ext) for ext in ['.txt', '.csv', '.tsv', '.h5', '.h5ad']):
            processed_data = read_scp_exp(unprocessed_data)
            success = True
        elif file_key.endswith('/'):
            # cellranger output 처리
            processed_data = read_10x_mtx(file_key)
            success = True
            
        if success:
            parts = file_key.split('/')
            information = {
                'study': parts[1],
                'type': parts[2],
                'file_name': parts[3],
                'data': processed_data
            }
            file_list.append(information)
        
    return file_list

def process_cluster(bucket_name, folder_prefix):
    response = list_s3_objects(bucket_name, folder_prefix)
        
    if response is None or 'Contents' not in response:
        # S3 객체가 없거나, 응답에 Contents가 없으면 빈 리스트 반환
        return []
        
    file_list = []
        
    for file in response['Contents']:
        file_key = file['Key']
        unprocessed_data = get_s3_objects(bucket_name, file_key)
        success = False
            
        # Check file type and process accordingly
        if any(file_key.endswith(ext) for ext in ['.txt', '.csv']):
            processed_data = read_scp_cluster(unprocessed_data)
            success = True
            
        if success:
            parts = file_key.split('/')
            information = {
                'study': parts[1],
                'type': parts[2],
                'file_name': parts[3],
                'data': processed_data
            }
            file_list.append(information)
        
    return file_list

def read_singlecellportal_from_s3(bucket_name, prefix):
    # Initialize file_list once and append results from both folders
    file_list = []
    file_list.extend(process_exp(bucket_name, prefix + 'expression/'))
    file_list.extend(process_cluster(bucket_name, prefix + 'cluster/'))
    
    return file_list

import tempfile

# 각 study 폴더 안에 있는 .h5ad 파일들을 리스트로 가져오는 함수
def get_h5ad_files_from_study(study_path):
    """
    특정 study 경로 안에 있는 모든 .h5ad 파일을 리스트로 가져옵니다.
    
    Parameters:
    - study_path: S3 버킷에서 study 폴더의 경로
    
    Returns:
    - h5ad_file_list: study 폴더 안의 모든 .h5ad 파일 경로 리스트
    """
    # study 폴더 안의 .h5ad 파일 리스트 초기화
    h5ad_file_list = []
    
    # S3 버킷의 해당 study 경로에서 파일들을 가져옴
    response_files = list_s3_objects('cellinsight-bucket', study_path, delimiter=True)
    
    print('response_files:', response_files)
    
    if response_files and 'Contents' in response_files:
        for obj in response_files['Contents']:
            file_key = obj['Key']
            print('file_key :', file_key)
            if file_key.endswith('.h5ad'):  # .h5ad 파일만 추가
                try:
                    # 임시 파일 생성
                    temp_file = tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False)
                
                    # S3에서 파일 다운로드
                    download_s3_objects('cellinsight-bucket', file_key, temp_file.name)
                
                    # 다운로드한 로컬 파일 경로를 리스트에 추가
                    h5ad_file_list.append(temp_file.name)
                
                except (OSError, IOError) as e:
                    # 오류 발생 시 로그를 남기고 해당 파일 건너뛰기
                    print(f"Error downloading or reading {file_key}: {e}")
                    # 임시 파일 삭제
                    temp_file.close()
                    os.remove(temp_file.name)
                    continue
    
    return h5ad_file_list