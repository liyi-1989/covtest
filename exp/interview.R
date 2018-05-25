# Complete the function below.

def  time_to_deliver(num_packages, delivery_sequence):
  n=num_packages
  major_grid=delivery_sequence[n-1][:4]
  count=0
  for i in range(n):
    if delivery_sequence[i][:4]==major_grid:
    count=count+1
  # here count=number of delivery
  # need to add number of move
  count=count+int(delivery_sequence[n-1][5:])-1
  for i in range(n-1):
    if (delivery_sequence[i][:4]!=major_grid) and (delivery_sequence[i+1][:4]==major_grid):
    major_num=int(delivery_sequence[i+1][5:])
  sub_num=int(delivery_sequence[i][5:])
  if (sub_num>=major_num):
    count=count+2*(sub_num-major_num)+1
  return(count)