import os
import pandas as pd

from django.utils.http import urlencode
from django.shortcuts import render, redirect
from django.http import HttpResponseRedirect

import preprocessing.umap as umap
from preprocessing import PanglaoDB_proc_python as pdl
from preprocessing import Single_Cell_Portal_proc as scp
from preprocessing.preproc import *
from processing.proc import *

import search.marker_search as marker_search
import search.gene_search as gene_search

import interaction.forms as forms

from . import functions

def only_render(request, html):
    return render(request, html)

def search(request):
    return render(request, 'search.html')

def preprocessing(request):
    if request.method == 'POST':
        form = forms.UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            file_format = form.cleaned_data['file_format']
            file_path = handle_uploaded_file(file)
            qc_form = forms.QCForm()
            return render(request, 'preprocessing.html', {
                'file_uploaded': True,
                'form': form,
                'qc_form': qc_form,
                'file_path': file_path,
                'file_format': file_format
            })
    else:
        form = forms.UploadFileForm()
        qc_form = forms.QCForm()
        return render(request, 'preprocessing.html', {
            'file_uploaded': False,
            'form': form,
            'qc_form': qc_form
        })

def handle_uploaded_file(f):
    media_path = 'media'
    if not os.path.exists(media_path):
        os.makedirs(media_path)
    file_path = os.path.join(media_path, f.name)
    with open(file_path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
    return file_path

def qc_process(request):
    if request.method == 'POST':
        qc_form = forms.QCForm(request.POST)
        file_path = request.POST.get('file_path')
        file_format = request.POST.get('file_format')

        if qc_form.is_valid():
            min_counts = qc_form.cleaned_data['min_counts']
            min_genes = qc_form.cleaned_data['min_genes']
            max_genes = qc_form.cleaned_data['max_genes']
            pct_counts_mt = qc_form.cleaned_data['pct_counts_mt']

            data = load_data(file_path, file_format)
            preprocessed_data = preprocess_data(data, min_counts, min_genes, max_genes, pct_counts_mt)
            violin_plot = draw_and_save_violin_plot(preprocessed_data)

            return render(request, 'preprocessing.html', {
                'violin_plot': violin_plot,
                'qc_form': qc_form,
                'show_mapcell_button': True,
                'file_uploaded': True,
                'file_path': file_path,
                'file_format': file_format,
                'min_counts': min_counts,
                'min_genes': min_genes,
                'max_genes': max_genes,
                'pct_counts_mt': pct_counts_mt
            })

    return HttpResponseRedirect('/preprocessing/')

def mapcell_process(request):
    if request.method == 'POST':
        file_path = request.POST.get('file_path')
        file_format = request.POST.get('file_format')
        min_counts = float(request.POST.get('min_counts'))
        min_genes = float(request.POST.get('min_genes'))
        max_genes = float(request.POST.get('max_genes'))
        pct_counts_mt = float(request.POST.get('pct_counts_mt'))

        # Normalization 및 UMAP Plot 생성

        data = load_data(file_path, file_format)
        preprocessed_data = preprocess_data(data, min_counts, min_genes, max_genes, pct_counts_mt)
        normalized_data = normalize_data(preprocessed_data)
        umap_plot = draw_and_save_umap_plot(normalized_data)

        # UMAP 플롯 페이지로 리디렉션
        umap_plot = umap_plot.replace('\\', '/')
        query_params = urlencode({'umap_plot': umap_plot})

        # GET 파라미터로 전달하여 UMAP 플롯 페이지로 리디렉션
        return redirect(f'/umap/?{query_params}')

    return HttpResponseRedirect('/preprocessing/')

def umap_view(request):
    user_umap_plot = request.GET.get('umap_plot')

    #Initialize
    s3_umap_plot, panglao_umap_plot = None, None
    csv_data, uns_data = None, None
    
    if request.method == 'POST':
        folder_name = request.POST.get('s3_file')
        delimiter = request.POST.get('delimiter', '\t')

        if folder_name:
            try:
                if 'PanglaoDB' in folder_name: 
                    #Select study of PanglaoDB

                    panglao_umap_plot, obs_data, uns_data = umap.fetch_and_process_file_Panglao(folder_name)
                    csv_data = pdl.import_additional_data(folder_name)
                else: 
                    #Select study of Singlecellportal

                    prefix_list = functions.load_cluster_file_list(folder_name)
                    cluster_file = prefix_list[0]
                    
                    if cluster_file:
                        s3_umap_plot = umap.fetch_and_process_file(cluster_file, delimiter)
                    else:
                        print("No cluster files found.")
            except Exception as e:
                print(f"Error fetching and processing file {folder_name}: {str(e)}")
                return render(request, 'umap.html', {
                    'user_umap_plot': user_umap_plot,
                    's3_umap_plot': None,
                    'panglao_umap_plot': None,
                    's3_files': functions.load_singlecellportal_folder_list(),
                    'csv_data': csv_data,
                    'uns_data': uns_data,
                    'error_message': f"Error processing file: {str(e)}"
                    
                })

    # S3 파일 목록 가져오기

    mapping = functions.load_singlecellportal_folder_list()
    PanglaoDB_files = functions.load_prefix_list('PanglaoDB/', delimiter = True)

    # singlecellportal 및 panglaoDB 파일들을 하나의 리스트로 합침
    combined_files = mapping + [(folder, folder) for folder in PanglaoDB_files]
    csv_data_list = csv_data.to_dict(orient='records') if csv_data is not None else []
    
    return render(request, 'umap.html', {
        'user_umap_plot': user_umap_plot,
        's3_umap_plot': s3_umap_plot,
        'panglao_umap_plot': panglao_umap_plot,
        's3_files': combined_files,  # 통합된 파일 리스트를 전달
        'csv_data': csv_data_list,
        'obs_data': obs_data if 'obs_data' in locals() else None,
        'uns_data': uns_data if 'uns_data' in locals() else None,
    })

    
from django.shortcuts import render
import search.marker_search as marker_search

def markersearch(request):
    selected_organ = request.GET.get('organ', 'All')
    selected_cell_type = request.GET.get('cell_type', 'All')

    # 필터링된 결과를 가져옴
    organ_cell_types, html_code = marker_search.fetch_markers(selected_organ, selected_cell_type)

    # 선택된 organ에 해당하는 cell types 리스트 가져오기
    cell_types = organ_cell_types.get(selected_organ, []) if selected_organ != 'All' else []

    # 템플릿으로 데이터 전달
    return render(request, 'markersearch.html', {
        'organ_cell_types_items': organ_cell_types.items(),
        'selected_organ': selected_organ,
        'selected_cell_type': selected_cell_type,
        'cell_types': cell_types,
        'html_code': html_code,
    })

from anndata import AnnData 

from migrations import models

def genesearch(request):
    if request.method == 'POST':
        form = forms.GeneSearchForm(request.POST)
        if form.is_valid():
            gene_name = form.cleaned_data['gene_name']
            
            # 초기화
            num_samples = 0
            num_clusters = 0
            study_list = []
            gene_existing_study_list = []

            # PanglaoDB 데이터셋 검색
            response = models.list_s3_objects('cellinsight-bucket', 'PanglaoDB/', delimiter=True)
            if response and 'CommonPrefixes' in response:
                study_list = [prefix['Prefix'] for prefix in response['CommonPrefixes']]
            
            for study in study_list:
                response_files = models.list_s3_objects('cellinsight-bucket', study, delimiter=True)
                study_file_list = []  # study_file_list 초기화
                if response_files and 'CommonPrefixes' in response_files:
                    study_file_list = [prefix['Prefix'] for prefix in response_files['CommonPrefixes']]
                    
                for file in study_file_list:
                    panglao_data = pdl.make_anndata(file)
                    if gene_search.check_gene_presence(panglao_data, gene_name):
                        num_samples += 1
                        # study 이름이 이미 리스트에 없으면 추가
                        if study not in gene_existing_study_list:
                            gene_existing_study_list.append(study)
                    
                    '''
                    # 클러스터 수를 얻어오고, 정수인지 확인한 후에 num_clusters에 더하기
                    cluster_count = gene_search.count_clusters_with_gene(adata, gene_name, cluster_key='leiden')
    
                    if isinstance(cluster_count, int):  # 클러스터 수가 정수일 경우
                        num_clusters += cluster_count
                    '''
            
            # SingleCellPortal 데이터셋 검색
            response = models.list_s3_objects('cellinsight-bucket', 'singlecellportal_anndata/', delimiter=True)
            if response and 'CommonPrefixes' in response:
                study_list = [prefix['Prefix'] for prefix in response['CommonPrefixes']]
                
            print('study_list:', study_list)
            
            for study in study_list:
                h5ad_files = models.get_h5ad_files_from_study(study)
                print('h5ad_files:', h5ad_files)
                for h5ad in h5ad_files:
                    try:
                        # anndata 파일 로드
                        if os.path.exists(h5ad) and h5ad.endswith('.h5ad'):
                            adata = sc.read_h5ad(h5ad)
                            
                            if gene_search.check_gene_presence(adata, gene_name):
                                num_samples += 1
                                # study 이름이 이미 리스트에 없으면 추가
                                if study not in gene_existing_study_list:
                                    gene_existing_study_list.append(study)
                    except (OSError, IOError) as e:
                        print(f"Error reading {h5ad}: {e}")
                        continue
                    
                    '''
                    # 클러스터 수를 얻어오고, 정수인지 확인한 후에 num_clusters에 더하기
                    cluster_count = gene_search.count_clusters_with_gene(adata, gene_name, cluster_key='leiden')
    
                    if isinstance(cluster_count, int):  # 클러스터 수가 정수일 경우
                        num_clusters += cluster_count
                    '''
            
            # 검색 결과를 템플릿에 전달하여 렌더링
            context = {
                'gene_name': gene_name,
                'num_samples': num_samples,
                'gene_existing_study_list': gene_existing_study_list
                #'num_clusters': num_clusters
            }
            return render(request, 'genesearch_result.html', context)
    else:
        form = forms.GeneSearchForm()

    return render(request, 'genesearch.html', {'form': form})