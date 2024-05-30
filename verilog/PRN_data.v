module primary_data_PRN_gen;
  function [0:54] R0;
    input [7:0] PRN_ID;
    begin
      case(PRN_ID)
        1: R0 = 55'o0061727026503255544;
        2: R0 = 55'o1660130752435362260;
        3: R0 = 55'o0676457016477551225;
        4: R0 = 55'o1763467705267605701;
        5: R0 = 55'o1614265052776007236;
        6: R0 = 55'o1446113457553463523;
        7: R0 = 55'o1467417471470124574;
        8: R0 = 55'o0022513456555401603;
        9: R0 = 55'o0004420115402210365;
        10: R0 = 55'o0072276243316574510;
        11: R0 = 55'o1632356715721616750;
        12: R0 = 55'o1670164755420300763;
        13: R0 = 55'o1752127524253360255;
        14: R0 = 55'o0262220014044243135;
        15: R0 = 55'o1476157654546440020;
        16: R0 = 55'o1567545246612304745;
        17: R0 = 55'o0341667641424721673;
        18: R0 = 55'o0627234635353763045;
        19: R0 = 55'o0422600144741165152;
        20: R0 = 55'o1661124176724621030;
        21: R0 = 55'o1225124173720602330;
        22: R0 = 55'o1271773065617322065;
        23: R0 = 55'o0611751161355750124;
        24: R0 = 55'o0121046615341766266;
        25: R0 = 55'o0337423707274604122;
        26: R0 = 55'o0246610305446052270;
        27: R0 = 55'o0427326063324033344;
        28: R0 = 55'o1127467544162733403;
        29: R0 = 55'o0772425336125565156;
        30: R0 = 55'o1652465113031101044;
        31: R0 = 55'o1737622607214524550;
        32: R0 = 55'o1621315362240732407;
        33: R0 = 55'o0171733204500613155;
        34: R0 = 55'o1462031354327077565;
        35: R0 = 55'o1141265411761074755;
        36: R0 = 55'o0665106277260231251;
        37: R0 = 55'o0573123144343776027;
        38: R0 = 55'o0222101406610314705;
        39: R0 = 55'o0140673225434336401;
        40: R0 = 55'o0624233245727625631;
        41: R0 = 55'o0224022145647544263;
        42: R0 = 55'o0222501602610354705;
        43: R0 = 55'o1370337660412244327;
        44: R0 = 55'o0563567347256715524;
        45: R0 = 55'o1407636661116077143;
        46: R0 = 55'o1137431557133151004;
        47: R0 = 55'o1113003456475500265;
        48: R0 = 55'o1746553632646152413;
        49: R0 = 55'o1465416631251321074;
        50: R0 = 55'o0130516430377202712;
        51: R0 = 55'o0762173527246302776;
        52: R0 = 55'o1606732407336425136;
        53: R0 = 55'o1131112010066741562;
        54: R0 = 55'o1107467740060732403;
        55: R0 = 55'o0755500241327076744;
        56: R0 = 55'o1443037764170374631;
        57: R0 = 55'o0243224434357700345;
        58: R0 = 55'o0445504023027564357;
        59: R0 = 55'o1211152271373271472;
        60: R0 = 55'o0256644102553071753;
        61: R0 = 55'o0733312314424771412;
        62: R0 = 55'o1636376400221406415;
        63: R0 = 55'o0574114621235461516;
        64: R0 = 55'o1710717574016037362;
      endcase
    end
  endfunction
  function [0:54] R1;
    input [7:0] PRN_ID;
    begin
      case(PRN_ID)
        1: R1 = 55'o0377627103341647600;
        2: R1 = 55'o0047555332635133703;
        3: R1 = 55'o0570574070736102152;
        4: R1 = 55'o0511013576745450615;
        5: R1 = 55'o1216243446624447775;
        6: R1 = 55'o0176452272675511054;
        7: R1 = 55'o0151055342317137706;
        8: R1 = 55'o1127720116046071664;
        9: R1 = 55'o0514407436155575524;
        10: R1 = 55'o0253070462740453542;
        11: R1 = 55'o0573371306324706336;
        12: R1 = 55'o1315135317732077306;
        13: R1 = 55'o1170303027726635012;
        14: R1 = 55'o1637171270537414673;
        15: R1 = 55'o0342370520251732111;
        16: R1 = 55'o0142423551056551362;
        17: R1 = 55'o0641261355426453710;
        18: R1 = 55'o0237176034757345266;
        19: R1 = 55'o1205663360515365064;
        20: R1 = 55'o0725000004121104102;
        21: R1 = 55'o0337367500320303262;
        22: R1 = 55'o1303374445022536530;
        23: R1 = 55'o1033071464007363115;
        24: R1 = 55'o0753124124237073577;
        25: R1 = 55'o0133522075443754772;
        26: R1 = 55'o1244212514312345145;
        27: R1 = 55'o1066056211234322164;
        28: R1 = 55'o0073115240113351010;
        29: R1 = 55'o1102260031574577224;
        30: R1 = 55'o1166703527236520553;
        31: R1 = 55'o0056062273631723177;
        32: R1 = 55'o0141517013160576212;
        33: R1 = 55'o1644007677312431616;
        34: R1 = 55'o0201757033615262622;
        35: R1 = 55'o0357610362675720200;
        36: R1 = 55'o1637504174727237065;
        37: R1 = 55'o1510345507743707753;
        38: R1 = 55'o0540160763721100120;
        39: R1 = 55'o0406415410457500342;
        40: R1 = 55'o0707515543554212732;
        41: R1 = 55'o0140216674314371011;
        42: R1 = 55'o0445414471314273300;
        43: R1 = 55'o0120121661750263177;
        44: R1 = 55'o0477301251340044262;
        45: R1 = 55'o1157040657040363676;
        46: R1 = 55'o1222265021477405004;
        47: R1 = 55'o0314661556545362364;
        48: R1 = 55'o0177320240371640542;
        49: R1 = 55'o0735517310345570340;
        50: R1 = 55'o1367565551220511432;
        51: R1 = 55'o1274167141162675644;
        52: R1 = 55'o1543641015130470077;
        53: R1 = 55'o0640733734534576460;
        54: R1 = 55'o0216312531021205434;
        55: R1 = 55'o0050232164401566177;
        56: R1 = 55'o0702636370401726111;
        57: R1 = 55'o1733537351460015703;
        58: R1 = 55'o1523265651140460620;
        59: R1 = 55'o0607703231502460135;
        60: R1 = 55'o1757246242710445777;
        61: R1 = 55'o0464412467237572274;
        62: R1 = 55'o1050617751566552643;
        63: R1 = 55'o1041606123021052264;
        64: R1 = 55'o1335441345250455042;
      endcase
    end
  endfunction
  function [0:4] C;
    input [7:0] PRN_ID;
    begin
      case(PRN_ID)
        1: C = 5'b10100;
        2: C = 5'b10100;
        3: C = 5'b00110;
        4: C = 5'b10100;
        5: C = 5'b10100;
        6: C = 5'b00110;
        7: C = 5'b10100;
        8: C = 5'b00110;
        9: C = 5'b00110;
        10: C = 5'b00110;
        11: C = 5'b10100;
        12: C = 5'b00110;
        13: C = 5'b10100;
        14: C = 5'b00110;
        15: C = 5'b00110;
        16: C = 5'b10100;
        17: C = 5'b00110;
        18: C = 5'b00110;
        19: C = 5'b00110;
        20: C = 5'b00110;
        21: C = 5'b10100;
        22: C = 5'b10100;
        23: C = 5'b10100;
        24: C = 5'b00110;
        25: C = 5'b10100;
        26: C = 5'b00110;
        27: C = 5'b00110;
        28: C = 5'b00110;
        29: C = 5'b00110;
        30: C = 5'b10100;
        31: C = 5'b10100;
        32: C = 5'b00110;
        33: C = 5'b10100;
        34: C = 5'b00110;
        35: C = 5'b00110;
        36: C = 5'b00110;
        37: C = 5'b10100;
        38: C = 5'b10100;
        39: C = 5'b01100;
        40: C = 5'b00110;
        41: C = 5'b00011;
        42: C = 5'b01100;
        43: C = 5'b10100;
        44: C = 5'b00110;
        45: C = 5'b10100;
        46: C = 5'b10100;
        47: C = 5'b00110;
        48: C = 5'b00110;
        49: C = 5'b00110;
        50: C = 5'b10100;
        51: C = 5'b10100;
        52: C = 5'b10100;
        53: C = 5'b00110;
        54: C = 5'b10100;
        55: C = 5'b00110;
        56: C = 5'b10100;
        57: C = 5'b00110;
        58: C = 5'b00110;
        59: C = 5'b10100;
        60: C = 5'b10010;
        61: C = 5'b10001;
        62: C = 5'b11000;
        63: C = 5'b00110;
        64: C = 5'b10100;
      endcase
    end
  endfunction
  function [0:10229] primary_data_PRN;  
    input [0:54] R0, R1;
    input [0:4] C;
    integer i;
    reg R0_fb, sigma_2A, sigma_2B, sigma_2C, sigma_2, R1A, R1B, R1_fb, C_fb;
    begin
      for(i=0; i<10230;i=i+1) begin
        R0_fb = R0[50]^R0[45]^R0[40]^R0[20]^R0[10]^R0[5]^R0[0];
        sigma_2A = (R0[50]^R0[45]^R0[40])&(R0[20]^R0[10]^R0[5]^R0[0]);
        sigma_2B = ((R0[50]^R0[45])&R0[40])^((R0[20]^R0[10])&(R0[5]^R0[0]));
        sigma_2C = (R0[50]&R0[45])^(R0[20]&R0[10])^(R0[5]&R0[0]);
        sigma_2 = sigma_2A^sigma_2B^sigma_2C;
        R1A = sigma_2 ^ (R0[40]^R0[35]^R0[30]^R0[25]^R0[15]^R0[0]);
        R1B = R1[50]^ R1[45]^R1[40]^ R1[20]^ R1[10]^ R1[5]^ R1[0];
        R1_fb = R1A ^ R1B;
        primary_data_PRN[i] = R1[0]^C[0];
        R0[0:53] =  R0[1:54]; 
        R0[54] = R0_fb;
        R1[0:53] = R1[1:54];
        R1[54] = R1_fb;
        C_fb = C[0];
        C[0:3] = C[1:4];
        C[4] = C_fb;
      end
      $display("First 24: %o",primary_data_PRN[0:23]);
      $display("Last 24: %o",primary_data_PRN[10206:10229]);
    end
  endfunction
  reg [7:0] Sat_num;
  reg[0:54] R0_init_val;
  reg[0:54] R1_init_val;
  reg[0:4] C_init_val;
  reg[0:10229] fin_PRN;
  initial begin
    Sat_num = 1;
    $display("PRN ID:%d",Sat_num);
    R0_init_val = R0(Sat_num);
    $display("INITIAL R0 Val: %o",R0_init_val);
    R1_init_val = R1(Sat_num);
    $display("INITIAL R1 Val: %o",R1_init_val);
    C_init_val = C(Sat_num);
    $display("INITIAL C Val: %b",C_init_val);
    fin_PRN = primary_data_PRN(R0_init_val, R1_init_val, C_init_val);
  end
endmodule
