function [Total_Number_of_Customers, Customers_Unserved] = update_customers_unserved_based_on_newArrivals(t, Nodes, ws, Total_Number_of_Customers, New_Arrivals_with_OD_t, Customers_Unserved)
New_Arrivals_with_OD_temp = New_Arrivals_with_OD_t{t, 1};

if size(New_Arrivals_with_OD_temp, 1) ~= 0
   No = (Total_Number_of_Customers + 1 : Total_Number_of_Customers + size(New_Arrivals_with_OD_temp, 1))';
   EnterTime = repmat(t, size(New_Arrivals_with_OD_temp, 1), 1);
   CTNode_No = [];
   Xcoor = [];
   Ycoor = [];
   D_x = [];
   D_y = [];
   CL = [];
   DrivingDis = [];
   BookingTime = [];
   for i = 1 : size(New_Arrivals_with_OD_temp, 1)
       CTNode_No = [CTNode_No; New_Arrivals_with_OD_temp(i, 1)];
       Xcoor = [Xcoor; New_Arrivals_with_OD_temp(i, 2)];
       Ycoor = [Ycoor; New_Arrivals_with_OD_temp(i, 3)];
       D_x = [D_x; New_Arrivals_with_OD_temp(i, 4)];
       D_y = [D_y; New_Arrivals_with_OD_temp(i, 5)];
       CL = [CL; Nodes(New_Arrivals_with_OD_temp(i, 1), 3)];
       DrivingDis = [DrivingDis; New_Arrivals_with_OD_temp(i, 7)];
       BookingTime = [BookingTime; New_Arrivals_with_OD_temp(i, 8)];
   end
   WaitingTime = zeros(size(New_Arrivals_with_OD_temp, 1), 1);
   WalkingSpeed = ws / 60 * ones(size(New_Arrivals_with_OD_temp, 1), 1); % walking speed is ws km/h = ws/60 km/min
   ReservedActive = zeros(size(New_Arrivals_with_OD_temp, 1), 1);
   ReservedPassive = zeros(size(New_Arrivals_with_OD_temp, 1), 1);
   Time2EV = NaN(size(New_Arrivals_with_OD_temp, 1), 1); % when a customer actively reserves an ev, he will walk to the ev to pick it up, this indicates the time when he arrives the ev
   Customers_Unserved_temp = [No, EnterTime, CTNode_No, Xcoor, Ycoor, D_x, D_y, CL, DrivingDis, BookingTime, WaitingTime, WalkingSpeed, ReservedActive, ReservedPassive, Time2EV];
else
   Customers_Unserved_temp = [];
end

Total_Number_of_Customers = Total_Number_of_Customers + size(New_Arrivals_with_OD_temp, 1);
Customers_Unserved{t, 1} = [Customers_Unserved{t, 1}; Customers_Unserved_temp];

if t >= 3
    Customers_Unserved{t - 2, 1} = [];
end