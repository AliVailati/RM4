function cf = coupons(coupon_dates, coupon_rate)
    %INPUTS: 
    %coupon_dates: a vector of coupon payment dates
    %coupon rate : a vector with the coupon rate for each period
    %OUTPUTS:
    %cf: a matrix with the coupon payment dates in the first column and the coupon payment in the other columns

    cf = zeros(length(coupon_dates), length(coupon_rate) + 1);  
    cf(:, 1) = coupon_dates;
    cf(:, 2:end) = ones(length(coupon_dates), 1)*coupon_rate;
    %Now we have to add one at the end of the matrix to all the coupons to consider the principal
    cf(end, 2:end) = cf(end, 2:end) + 1;
end