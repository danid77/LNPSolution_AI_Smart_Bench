// datatables-config.js
export function initializeDataTable(tableSelector, options = {}) {
    // 기본 옵션 정의
    const defaultOptions = {
        responsive: true,
        paging: true, // 페이지네이션 활성화
        lengthChange: true, // 표시 행 개수 선택 활성화
        pageLength: 10, // 기본 표시 행 개수
        lengthMenu: [10, 50, 100, 200, 300], // 선택 가능한 표시 행 개수
        pagingType: 'full_numbers', // 페이지네이션 버튼 유형 (숫자형)
        dom: '<"row"<"col-md-6 d-flex align-items-center"l><"col-md-6"f>>tp', // 검색창과 표시 개수 선택 정렬
        language: {
            search: "검색:",
            lengthMenu: "_MENU_ 항목 표시",
            paginate: {
                first: "처음",
                last: "마지막",
                next: "다음",
                previous: "이전"
            },
            info: "총 _TOTAL_개의 데이터 중 _START_ ~ _END_",
        },
        columnDefs: [
            { targets: "_all", className: "dt-center" } // 텍스트 가운데 정렬
        ]
    };

    // 사용자 지정 옵션과 병합
    const finalOptions = { ...defaultOptions, ...options };

    // DataTables 초기화
    const table = $(tableSelector).DataTable(finalOptions);

    // 추가 스타일링
    $('.dataTables_filter input').addClass('form-control').css({
        borderRadius: '5px',
        marginLeft: '0.5rem',
        width: '300px', // 검색창 가로 크기 증가
        maxWidth: '100%', // 반응형 제한
    });

    $('.dataTables_length select').addClass('form-select').css({
        width: 'auto', // 드롭다운 길이 자동 조정
        padding: '5px', // 패딩 추가
        minWidth: '150px' // 최소 가로 길이 설정
    });

    // 페이지네이션 가운데 정렬
    $('.dataTables_paginate').addClass('d-flex justify-content-center mt-3');
    
    // 테이블 데이터 수 표시
    const tableInfo = table.page.info();
    const totalRows = tableInfo.recordsTotal; // 전체 데이터 수
    const message = `현재 테이블에는 총 ${totalRows}개의 데이터가 있습니다.`;
    console.log(message);

    // 화면에 표시 (옵션)
    $(tableSelector).before(`<div class="table-info mb-3 text-center">${message}</div>`);

    return table; // 초기화된 DataTable 객체 반환
}

export function downloadTables(fixedHeaders, dataFrame, accession, jobName, saveFileName, options = {}) {

    if (dataFrame && dataFrame.length > 0) {
        // 고정된 컬럼 순서 : 외부에서 가져올거임
        // const fixedHeaders = ['Label', 'Antisense', 'Region1', 'Region2', 'Exon', 'CDS', 'GC', 'Score'];
    
        // 데이터의 모든 컬럼 가져오기
        const allHeaders = Object.keys(dataFrame[0]);
    
        // 고정된 컬럼을 제외한 나머지 컬럼 추가
        const additionalHeaders = allHeaders.filter(header => !fixedHeaders.includes(header));
    
        // 최종 헤더 배열 생성
        const finalHeaders = [...fixedHeaders, ...additionalHeaders];

        if (!Array.isArray(dataFrame) || dataFrame.length === 0) {
            console.error('Invalid or empty dataFrame');
            alert('No valid data available to download.');
            return;
        }

        const csvContent = [
            finalHeaders.join(','), // Add headers
            ...dataFrame.map(row => finalHeaders.map(header => `"${row[header] || ''}"`).join(',')) // Add rows
        ].join('\n');

        // Create Blob for the CSV file
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const url = URL.createObjectURL(blob);

        const link = document.createElement('a');
        link.href = url;
        link.download = `${accession}_${jobName}_${saveFileName}.csv`; // File name
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }else {
        alert('No data available for download.');
    }
}