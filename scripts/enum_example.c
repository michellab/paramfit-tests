//Test on enum

enum day { sun, mon, tue, wed, thu, fri, sat};

typedef enum day   day;

day find_next_day(day d)
{
  day next_day;

  return ((day)(((int) d+ + 1) % 7));
}


int main ()
{
  day  which;
  which = mon;
  find_next_day(which);
}
