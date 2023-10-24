function binaryVector = convertToBinary(vector)
    binaryVector = zeros(size(vector));  % Khởi tạo vector kết quả với giá trị 0 ban đầu
    binaryVector(vector >= 0.5) = 1;  % Thiết lập các phần tử lớn hơn hoặc bằng 0.5 thành 1
end