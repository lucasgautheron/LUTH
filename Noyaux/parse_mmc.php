<?php
$lines = file('mmc01.txt');

$A = $Z = $Q = 0;

$fp = fopen("single_rates", "w+");
foreach($lines as $line)
{
    $values = array();
    if(preg_match("#z=\s*(.*) n=\s*(.*) a=\s*(.*) Q=\s*(.*)#", $line, $values))
    {
        $Z = $values[1];
        $A = $values[3];
        $Q = $values[4];
    }
    else if(!preg_match("#[a-zA-Z]#", $line))
    {
        $str = "\t$A\t$Z\t$Q\t$line";
        echo $str;
        fwrite($fp, $str);
    }
}
fclose($fp);
?>
