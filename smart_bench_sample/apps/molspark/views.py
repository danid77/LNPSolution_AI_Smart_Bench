from django.shortcuts import render
from django.http import JsonResponse
from django.http import HttpResponse
import json, os

from apps.molspark import services as sv
# Create your views here.

def molsparkInputPage(request):
    return render(request, 'molspark/input.html')

def molsparkOutputPage(request):
    return render(request, 'molspark/output.html')


def molsparkCal(request):
    try:
        # JSON 데이터 파싱
        data = json.loads(request.body)
        smiles = data.get("smiles")
        model_type = data.get("model_type")
        strategy = data.get("strategy")
        temperature = data.get("temperature")
        num_samples = data.get("num_samples")
        
        # print(smiles, model_type, strategy, temperature, num_samples)
        molspark_result = sv.molsparkProcess(smiles, model_type, strategy, temperature, num_samples)
        print(molspark_result)
        return JsonResponse({"status": "success", "molspark_result": molspark_result})
    except Exception as e:
        return JsonResponse({"status": "error", "message": str(e)})


def molsparkTable(request):
    """ JSON 데이터를 반환하는 API 뷰 """
    file_path = request.GET.get('result', '')

    # 파일 경로 로그 출력 (디버깅용)
    print(f"📂 Received file path: {file_path}")

    if not os.path.isfile(file_path):
        print("❌ 파일이 존재하지 않습니다!")
        return JsonResponse({"status": "error", "message": "File not found"}, status=404)

    try:
        dict_data = sv.csv_to_dict(file_path)  # CSV를 딕셔너리 형태로 변환

        if dict_data is None:
            print("❌ CSV 파일을 변환하는 중 오류 발생!")
            return JsonResponse({"status": "error", "message": "CSV parsing error"}, status=500)

        print(f"✅ CSV 변환 성공, 데이터 개수: {len(dict_data)} 개")
        return JsonResponse({"status": "success", "data": dict_data})  # JSON 응답 반환
    except Exception as e:
        print(f"❌ 데이터 처리 중 예외 발생: {e}")
        return JsonResponse({"status": "error", "message": str(e)}, status=500)

def test(request):
    return render(request, 'molspark.html')