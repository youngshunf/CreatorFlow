/**
 * China Administrative Regions Data
 *
 * Hierarchical data structure for provinces, cities, and districts.
 */

export interface Region {
  code: string
  name: string
}

export interface City extends Region {
  districts: Region[]
}

export interface Province extends Region {
  cities: City[]
}

/**
 * China regions data: provinces, cities, and districts
 */
export const chinaRegions: Province[] = [
  {
    code: '110000',
    name: '北京市',
    cities: [
      {
        code: '110100',
        name: '北京市',
        districts: [
          { code: '110101', name: '东城区' },
          { code: '110102', name: '西城区' },
          { code: '110105', name: '朝阳区' },
          { code: '110106', name: '丰台区' },
          { code: '110107', name: '石景山区' },
          { code: '110108', name: '海淀区' },
          { code: '110109', name: '门头沟区' },
          { code: '110111', name: '房山区' },
          { code: '110112', name: '通州区' },
          { code: '110113', name: '顺义区' },
          { code: '110114', name: '昌平区' },
          { code: '110115', name: '大兴区' },
          { code: '110116', name: '怀柔区' },
          { code: '110117', name: '平谷区' },
          { code: '110118', name: '密云区' },
          { code: '110119', name: '延庆区' },
        ],
      },
    ],
  },
  {
    code: '120000',
    name: '天津市',
    cities: [
      {
        code: '120100',
        name: '天津市',
        districts: [
          { code: '120101', name: '和平区' },
          { code: '120102', name: '河东区' },
          { code: '120103', name: '河西区' },
          { code: '120104', name: '南开区' },
          { code: '120105', name: '河北区' },
          { code: '120106', name: '红桥区' },
          { code: '120110', name: '东丽区' },
          { code: '120111', name: '西青区' },
          { code: '120112', name: '津南区' },
          { code: '120113', name: '北辰区' },
          { code: '120114', name: '武清区' },
          { code: '120115', name: '宝坻区' },
          { code: '120116', name: '滨海新区' },
        ],
      },
    ],
  },
  {
    code: '310000',
    name: '上海市',
    cities: [
      {
        code: '310100',
        name: '上海市',
        districts: [
          { code: '310101', name: '黄浦区' },
          { code: '310104', name: '徐汇区' },
          { code: '310105', name: '长宁区' },
          { code: '310106', name: '静安区' },
          { code: '310107', name: '普陀区' },
          { code: '310109', name: '虹口区' },
          { code: '310110', name: '杨浦区' },
          { code: '310112', name: '闵行区' },
          { code: '310113', name: '宝山区' },
          { code: '310114', name: '嘉定区' },
          { code: '310115', name: '浦东新区' },
          { code: '310116', name: '金山区' },
          { code: '310117', name: '松江区' },
          { code: '310118', name: '青浦区' },
          { code: '310120', name: '奉贤区' },
          { code: '310151', name: '崇明区' },
        ],
      },
    ],
  },
  {
    code: '500000',
    name: '重庆市',
    cities: [
      {
        code: '500100',
        name: '重庆市',
        districts: [
          { code: '500101', name: '万州区' },
          { code: '500102', name: '涪陵区' },
          { code: '500103', name: '渝中区' },
          { code: '500104', name: '大渡口区' },
          { code: '500105', name: '江北区' },
          { code: '500106', name: '沙坪坝区' },
          { code: '500107', name: '九龙坡区' },
          { code: '500108', name: '南岸区' },
          { code: '500109', name: '北碚区' },
          { code: '500110', name: '綦江区' },
          { code: '500111', name: '大足区' },
          { code: '500112', name: '渝北区' },
          { code: '500113', name: '巴南区' },
        ],
      },
    ],
  },
  {
    code: '440000',
    name: '广东省',
    cities: [
      {
        code: '440100',
        name: '广州市',
        districts: [
          { code: '440103', name: '荔湾区' },
          { code: '440104', name: '越秀区' },
          { code: '440105', name: '海珠区' },
          { code: '440106', name: '天河区' },
          { code: '440111', name: '白云区' },
          { code: '440112', name: '黄埔区' },
          { code: '440113', name: '番禺区' },
          { code: '440114', name: '花都区' },
          { code: '440115', name: '南沙区' },
          { code: '440117', name: '从化区' },
          { code: '440118', name: '增城区' },
        ],
      },
      {
        code: '440300',
        name: '深圳市',
        districts: [
          { code: '440303', name: '罗湖区' },
          { code: '440304', name: '福田区' },
          { code: '440305', name: '南山区' },
          { code: '440306', name: '宝安区' },
          { code: '440307', name: '龙岗区' },
          { code: '440308', name: '盐田区' },
          { code: '440309', name: '龙华区' },
          { code: '440310', name: '坪山区' },
          { code: '440311', name: '光明区' },
        ],
      },
      {
        code: '440400',
        name: '珠海市',
        districts: [
          { code: '440402', name: '香洲区' },
          { code: '440403', name: '斗门区' },
          { code: '440404', name: '金湾区' },
        ],
      },
      {
        code: '440600',
        name: '佛山市',
        districts: [
          { code: '440604', name: '禅城区' },
          { code: '440605', name: '南海区' },
          { code: '440606', name: '顺德区' },
          { code: '440607', name: '三水区' },
          { code: '440608', name: '高明区' },
        ],
      },
      {
        code: '441300',
        name: '惠州市',
        districts: [
          { code: '441302', name: '惠城区' },
          { code: '441303', name: '惠阳区' },
          { code: '441322', name: '博罗县' },
          { code: '441323', name: '惠东县' },
          { code: '441324', name: '龙门县' },
        ],
      },
      {
        code: '441900',
        name: '东莞市',
        districts: [
          { code: '441900', name: '东莞市' },
        ],
      },
    ],
  },
  {
    code: '330000',
    name: '浙江省',
    cities: [
      {
        code: '330100',
        name: '杭州市',
        districts: [
          { code: '330102', name: '上城区' },
          { code: '330105', name: '拱墅区' },
          { code: '330106', name: '西湖区' },
          { code: '330108', name: '滨江区' },
          { code: '330109', name: '萧山区' },
          { code: '330110', name: '余杭区' },
          { code: '330111', name: '富阳区' },
          { code: '330112', name: '临安区' },
          { code: '330113', name: '临平区' },
          { code: '330114', name: '钱塘区' },
        ],
      },
      {
        code: '330200',
        name: '宁波市',
        districts: [
          { code: '330203', name: '海曙区' },
          { code: '330205', name: '江北区' },
          { code: '330206', name: '北仑区' },
          { code: '330211', name: '镇海区' },
          { code: '330212', name: '鄞州区' },
          { code: '330213', name: '奉化区' },
        ],
      },
      {
        code: '330300',
        name: '温州市',
        districts: [
          { code: '330302', name: '鹿城区' },
          { code: '330303', name: '龙湾区' },
          { code: '330304', name: '瓯海区' },
          { code: '330305', name: '洞头区' },
        ],
      },
    ],
  },
  {
    code: '320000',
    name: '江苏省',
    cities: [
      {
        code: '320100',
        name: '南京市',
        districts: [
          { code: '320102', name: '玄武区' },
          { code: '320104', name: '秦淮区' },
          { code: '320105', name: '建邺区' },
          { code: '320106', name: '鼓楼区' },
          { code: '320111', name: '浦口区' },
          { code: '320113', name: '栖霞区' },
          { code: '320114', name: '雨花台区' },
          { code: '320115', name: '江宁区' },
        ],
      },
      {
        code: '320500',
        name: '苏州市',
        districts: [
          { code: '320505', name: '虎丘区' },
          { code: '320506', name: '吴中区' },
          { code: '320507', name: '相城区' },
          { code: '320508', name: '姑苏区' },
          { code: '320509', name: '吴江区' },
          { code: '320581', name: '常熟市' },
          { code: '320582', name: '张家港市' },
          { code: '320583', name: '昆山市' },
          { code: '320585', name: '太仓市' },
        ],
      },
      {
        code: '320200',
        name: '无锡市',
        districts: [
          { code: '320205', name: '锡山区' },
          { code: '320206', name: '惠山区' },
          { code: '320211', name: '滨湖区' },
          { code: '320213', name: '梁溪区' },
          { code: '320214', name: '新吴区' },
          { code: '320281', name: '江阴市' },
          { code: '320282', name: '宜兴市' },
        ],
      },
    ],
  },
  {
    code: '510000',
    name: '四川省',
    cities: [
      {
        code: '510100',
        name: '成都市',
        districts: [
          { code: '510104', name: '锦江区' },
          { code: '510105', name: '青羊区' },
          { code: '510106', name: '金牛区' },
          { code: '510107', name: '武侯区' },
          { code: '510108', name: '成华区' },
          { code: '510112', name: '龙泉驿区' },
          { code: '510113', name: '青白江区' },
          { code: '510114', name: '新都区' },
          { code: '510115', name: '温江区' },
          { code: '510116', name: '双流区' },
          { code: '510117', name: '郫都区' },
          { code: '510118', name: '新津区' },
        ],
      },
    ],
  },
  {
    code: '420000',
    name: '湖北省',
    cities: [
      {
        code: '420100',
        name: '武汉市',
        districts: [
          { code: '420102', name: '江岸区' },
          { code: '420103', name: '江汉区' },
          { code: '420104', name: '硚口区' },
          { code: '420105', name: '汉阳区' },
          { code: '420106', name: '武昌区' },
          { code: '420107', name: '青山区' },
          { code: '420111', name: '洪山区' },
          { code: '420112', name: '东西湖区' },
          { code: '420113', name: '汉南区' },
          { code: '420114', name: '蔡甸区' },
          { code: '420115', name: '江夏区' },
          { code: '420116', name: '黄陂区' },
          { code: '420117', name: '新洲区' },
        ],
      },
    ],
  },
  {
    code: '430000',
    name: '湖南省',
    cities: [
      {
        code: '430100',
        name: '长沙市',
        districts: [
          { code: '430102', name: '芙蓉区' },
          { code: '430103', name: '天心区' },
          { code: '430104', name: '岳麓区' },
          { code: '430105', name: '开福区' },
          { code: '430111', name: '雨花区' },
          { code: '430112', name: '望城区' },
          { code: '430121', name: '长沙县' },
        ],
      },
    ],
  },
  {
    code: '350000',
    name: '福建省',
    cities: [
      {
        code: '350100',
        name: '福州市',
        districts: [
          { code: '350102', name: '鼓楼区' },
          { code: '350103', name: '台江区' },
          { code: '350104', name: '仓山区' },
          { code: '350105', name: '马尾区' },
          { code: '350111', name: '晋安区' },
          { code: '350112', name: '长乐区' },
        ],
      },
      {
        code: '350200',
        name: '厦门市',
        districts: [
          { code: '350203', name: '思明区' },
          { code: '350205', name: '海沧区' },
          { code: '350206', name: '湖里区' },
          { code: '350211', name: '集美区' },
          { code: '350212', name: '同安区' },
          { code: '350213', name: '翔安区' },
        ],
      },
    ],
  },
  {
    code: '370000',
    name: '山东省',
    cities: [
      {
        code: '370100',
        name: '济南市',
        districts: [
          { code: '370102', name: '历下区' },
          { code: '370103', name: '市中区' },
          { code: '370104', name: '槐荫区' },
          { code: '370105', name: '天桥区' },
          { code: '370112', name: '历城区' },
          { code: '370113', name: '长清区' },
          { code: '370114', name: '章丘区' },
          { code: '370115', name: '济阳区' },
          { code: '370116', name: '莱芜区' },
          { code: '370117', name: '钢城区' },
        ],
      },
      {
        code: '370200',
        name: '青岛市',
        districts: [
          { code: '370202', name: '市南区' },
          { code: '370203', name: '市北区' },
          { code: '370211', name: '黄岛区' },
          { code: '370212', name: '崂山区' },
          { code: '370213', name: '李沧区' },
          { code: '370214', name: '城阳区' },
          { code: '370215', name: '即墨区' },
        ],
      },
    ],
  },
  {
    code: '410000',
    name: '河南省',
    cities: [
      {
        code: '410100',
        name: '郑州市',
        districts: [
          { code: '410102', name: '中原区' },
          { code: '410103', name: '二七区' },
          { code: '410104', name: '管城回族区' },
          { code: '410105', name: '金水区' },
          { code: '410106', name: '上街区' },
          { code: '410108', name: '惠济区' },
        ],
      },
    ],
  },
  {
    code: '610000',
    name: '陕西省',
    cities: [
      {
        code: '610100',
        name: '西安市',
        districts: [
          { code: '610102', name: '新城区' },
          { code: '610103', name: '碑林区' },
          { code: '610104', name: '莲湖区' },
          { code: '610111', name: '灞桥区' },
          { code: '610112', name: '未央区' },
          { code: '610113', name: '雁塔区' },
          { code: '610114', name: '阎良区' },
          { code: '610115', name: '临潼区' },
          { code: '610116', name: '长安区' },
          { code: '610117', name: '高陵区' },
          { code: '610118', name: '鄠邑区' },
        ],
      },
    ],
  },
  {
    code: '130000',
    name: '河北省',
    cities: [
      {
        code: '130100',
        name: '石家庄市',
        districts: [
          { code: '130102', name: '长安区' },
          { code: '130104', name: '桥西区' },
          { code: '130105', name: '新华区' },
          { code: '130107', name: '井陉矿区' },
          { code: '130108', name: '裕华区' },
          { code: '130109', name: '藁城区' },
          { code: '130110', name: '鹿泉区' },
          { code: '130111', name: '栾城区' },
        ],
      },
    ],
  },
  {
    code: '210000',
    name: '辽宁省',
    cities: [
      {
        code: '210100',
        name: '沈阳市',
        districts: [
          { code: '210102', name: '和平区' },
          { code: '210103', name: '沈河区' },
          { code: '210104', name: '大东区' },
          { code: '210105', name: '皇姑区' },
          { code: '210106', name: '铁西区' },
          { code: '210111', name: '苏家屯区' },
          { code: '210112', name: '浑南区' },
          { code: '210113', name: '沈北新区' },
          { code: '210114', name: '于洪区' },
        ],
      },
      {
        code: '210200',
        name: '大连市',
        districts: [
          { code: '210202', name: '中山区' },
          { code: '210203', name: '西岗区' },
          { code: '210204', name: '沙河口区' },
          { code: '210211', name: '甘井子区' },
          { code: '210212', name: '旅顺口区' },
          { code: '210213', name: '金州区' },
          { code: '210214', name: '普兰店区' },
        ],
      },
    ],
  },
  {
    code: '340000',
    name: '安徽省',
    cities: [
      {
        code: '340100',
        name: '合肥市',
        districts: [
          { code: '340102', name: '瑶海区' },
          { code: '340103', name: '庐阳区' },
          { code: '340104', name: '蜀山区' },
          { code: '340111', name: '包河区' },
          { code: '340121', name: '长丰县' },
          { code: '340122', name: '肥东县' },
          { code: '340123', name: '肥西县' },
          { code: '340124', name: '庐江县' },
        ],
      },
    ],
  },
  {
    code: '360000',
    name: '江西省',
    cities: [
      {
        code: '360100',
        name: '南昌市',
        districts: [
          { code: '360102', name: '东湖区' },
          { code: '360103', name: '西湖区' },
          { code: '360104', name: '青云谱区' },
          { code: '360111', name: '青山湖区' },
          { code: '360112', name: '新建区' },
          { code: '360113', name: '红谷滩区' },
        ],
      },
    ],
  },
  {
    code: '450000',
    name: '广西壮族自治区',
    cities: [
      {
        code: '450100',
        name: '南宁市',
        districts: [
          { code: '450102', name: '兴宁区' },
          { code: '450103', name: '青秀区' },
          { code: '450105', name: '江南区' },
          { code: '450107', name: '西乡塘区' },
          { code: '450108', name: '良庆区' },
          { code: '450109', name: '邕宁区' },
          { code: '450110', name: '武鸣区' },
        ],
      },
    ],
  },
  {
    code: '530000',
    name: '云南省',
    cities: [
      {
        code: '530100',
        name: '昆明市',
        districts: [
          { code: '530102', name: '五华区' },
          { code: '530103', name: '盘龙区' },
          { code: '530111', name: '官渡区' },
          { code: '530112', name: '西山区' },
          { code: '530113', name: '东川区' },
          { code: '530114', name: '呈贡区' },
          { code: '530115', name: '晋宁区' },
        ],
      },
    ],
  },
  {
    code: '520000',
    name: '贵州省',
    cities: [
      {
        code: '520100',
        name: '贵阳市',
        districts: [
          { code: '520102', name: '南明区' },
          { code: '520103', name: '云岩区' },
          { code: '520111', name: '花溪区' },
          { code: '520112', name: '乌当区' },
          { code: '520113', name: '白云区' },
          { code: '520115', name: '观山湖区' },
        ],
      },
    ],
  },
  {
    code: '460000',
    name: '海南省',
    cities: [
      {
        code: '460100',
        name: '海口市',
        districts: [
          { code: '460105', name: '秀英区' },
          { code: '460106', name: '龙华区' },
          { code: '460107', name: '琼山区' },
          { code: '460108', name: '美兰区' },
        ],
      },
      {
        code: '460200',
        name: '三亚市',
        districts: [
          { code: '460202', name: '海棠区' },
          { code: '460203', name: '吉阳区' },
          { code: '460204', name: '天涯区' },
          { code: '460205', name: '崖州区' },
        ],
      },
    ],
  },
]

/**
 * Get all provinces
 */
export function getProvinces(): Region[] {
  return chinaRegions.map(({ code, name }) => ({ code, name }))
}

/**
 * Get cities for a province
 */
export function getCitiesByProvince(provinceName: string): Region[] {
  const province = chinaRegions.find((p) => p.name === provinceName)
  if (!province) return []
  return province.cities.map(({ code, name }) => ({ code, name }))
}

/**
 * Get districts for a city
 */
export function getDistrictsByCity(provinceName: string, cityName: string): Region[] {
  const province = chinaRegions.find((p) => p.name === provinceName)
  if (!province) return []
  const city = province.cities.find((c) => c.name === cityName)
  if (!city) return []
  return city.districts
}
