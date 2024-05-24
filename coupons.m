function cf = coupons(coupon_dates, coupon_rate)
    cf = zeros(length(coupon_dates), 2);
    cf(:, 1) = coupon_dates;
    cf(:, 2) = coupon_rate;
    cf(end, end) = 1 + coupon_rate;
end