"""Hiroshima University parking SP: reproducible illustrative respondents, not observations.
Run: python 01_generate_data.py
Inputs: design.csv, simulation_settings.json. Standard library only.
"""
import csv, json, math, random, pathlib, collections

BASE = pathlib.Path(__file__).resolve().parent
CONFIG = json.loads((BASE/'simulation_settings.json').read_text(encoding='utf-8'))
RNG = random.Random(CONFIG['seed'])

def csv_read(name):
    with (BASE/name).open('r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def csv_write(name, rows, fields=None):
    rows=list(rows)
    fields=fields or list(rows[0])
    with (BASE/name).open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def expit(z):
    if z >= 0:
        return 1/(1+math.exp(-z))
    e=math.exp(z)
    return e/(1+e)

def clip(v, lo, hi): return min(hi, max(lo, v))
def sample_weighted(labels, weights): return RNG.choices(labels, weights=weights, k=1)[0]

# Explicit exact margins across each information group (n=200 per group).
# Student/parking-user=64, student/nonuser=96, faculty/staff-user=16, staff/nonuser=24.
info_map={'I1':'追加情報なし','I2':'ゲート更新','I3':'カーボンニュートラル','I4':'他大学料金','I5':'駐車場環境整備'}
people=[]
for group_i in range(1,6):
    ig=f'I{group_i}'
    combos = [('学生',1)]*64+[('学生',0)]*96+[('教職員',1)]*16+[('教職員',0)]*24
    RNG.shuffle(combos)
    # 67/67/66 rotated by information group -> balanced 15 versions.
    base=[67,67,66]; rotated=base[group_i-1:]+base[:group_i-1] if group_i<=3 else None
    # rotate cyclically (0,1,2,0,1) instead of slicing past bounds
    shift=(group_i-1)%3
    block_sizes=base[shift:]+base[:shift]
    blocks=[b for b,n in enumerate(block_sizes,1) for _ in range(n)]
    RNG.shuffle(blocks)
    for (segment,user),block in zip(combos,blocks):
        role=sample_weighted(['学部生','大学院生'],[0.86,0.14]) if segment=='学生' else sample_weighted(['教員','職員'],[0.34,0.66])
        if user:
            days=sample_weighted([1,2,3,4,5,6,7],[5,10,20,25,33,5,2])
            commute='自動車（自分で運転）' if RNG.random()<0.90 else sample_weighted(['バス','自転車','鉄道＋バス等'],[2,3,1])
            car_access=1
        else:
            days=0
            commute=sample_weighted(['自転車','徒歩','バス','鉄道＋バス等','自動車（他者が運転）','バイク'],[32,30,18,12,5,3])
            car_access=int(RNG.random() < (0.55 if segment=='教職員' else 0.25))
        commute_freq=sample_weighted(['週5日以上','週3～4日','週1～2日','月1～3日'],[65,20,12,3])
        # Illustrative district categories, NOT real or geocoded addresses.
        departure='架空設定：'+sample_weighted(['東広島市内','広島市内','呉市周辺','その他地域'], [77,13,7,3])
        # Baseline car use is the agreed user=drives at least one day/week.
        fee_awareness=sample_weighted(['知っており支払う','知っているが払わない','知らなかった'],[92,6,2]) if user else sample_weighted(['知っており支払う','知っているが払わない','知らなかった'],[2,48,50])
        fee_perception=sample_weighted(['非常に安い','やや安い','適切','やや高い','非常に高い','判断できない'],[14,35,35,10,4,2]) if user else sample_weighted(['非常に安い','やや安い','適切','やや高い','非常に高い','判断できない'],[9,26,27,12,5,21])
        knowledge=sample_weighted(['知っていた','一部知っていた','知らなかった'],[14,34,52]) if group_i>1 else ''
        reason_all=['所要時間','公共交通便数','送迎・介護','障害・健康上の事情','費用','勤務・授業時間','天候','駐車しやすさ','その他']
        reason_probs=[.78,.60,.08,.04,.45,.42,.23,.23,.06] if user else [.22,.24,.04,.03,.10,.12,.14,.07,.03]
        reasons={f'q6_{i+1}':int(RNG.random()<p) for i,p in enumerate(reason_probs)}
        if not any(reasons.values()): reasons['q6_1']=1  # mandatory in supplied form
        q9_opts={f'q9_{i+1}':int(RNG.random()<p) for i,p in enumerate([0.45,0.66,0.18,0.04])}
        if not any(q9_opts.values()): q9_opts['q9_5']=1
        else:q9_opts['q9_5']=0
        q10={f'q10_{i+1}': int(clip(round(RNG.gauss(3.1 if user else 2.65,0.95)),1,5)) for i in range(7)}
        e1_probs=[.50,.20,.15,.10,.30,.05] if user else [.32,.17,.12,.10,.26,.04]
        e1={f'e1_{i+1}':int(RNG.random()<p) for i,p in enumerate(e1_probs)}
        e1['e1_7']=int(not any(e1.values()))
        # Individual heterogeneity: shared across 6 tasks, retained in simulation only (excluded from observed data).
        ppl={
            'respondent_id':'R'+str(len(people)+1).zfill(4),'survey_version':f'{ig}-{block}',
            'info_group':ig,'info_label':info_map[ig],'block':block,
            'segment':segment,'role':role,'parking_user':user,'qs1_consent':1,'qs2_frequency':commute_freq,
            'q1_role':role,'q2_main_mode':commute,'q3_car_access':car_access,
            'q4_car_days_week':days,'q5_departure_text':departure,
            'q7_fee_awareness':fee_awareness,'q8_current_fee_view':fee_perception,
            'q11_prior_info':knowledge,'e2_free_text':'',
            '_attitude':RNG.gauss(0,0.36), '_price_sensitivity':clip(RNG.gauss(1,0.20),0.45,1.60),
            '_parking_propensity':RNG.gauss(0,0.43),
            **reasons,**q9_opts,**q10,**e1
        }
        people.append(ppl)
RNG.shuffle(people)
# Reassign IDs to shuffled respondents, to avoid ordering carrying treatment labels.
for idx,p in enumerate(people,1):p['respondent_id']='R'+str(idx).zfill(4)

# design.csv created from the user's 18-task table in the 2026-Nov memo.
design=csv_read('design.csv')
assert len(design)==54
by_block=collections.defaultdict(lambda: collections.defaultdict(list))
for d in design:
    d['task_id']=int(d['task_id']);d['block']=int(d['block'])
    for k in ['price_yen','gate_code','safety_code','parking_service_code','bus_discount_code','bus_frequency_code','is_reform']:
        d[k]=int(d[k])
    by_block[d['block']][d['task_id']].append(d)
for b in (1,2,3):
    assert len(by_block[b])==6
    for task,options in by_block[b].items():
        assert [opt['alternative'] for opt in options]==['A','B','SQ']

cf=CONFIG['choice_coefficients']
rf=CONFIG['support_parameters']
uf=CONFIG['parking_intention_parameters']
choice_rows=[]
task_rows=[]
for p in people:
    user=p['parking_user']; inf=p['info_group']; rid=p['respondent_id']
    for task,alts in by_block[p['block']].items():
        latent=[]
        for d in alts:
            reform=d['is_reform'];
            price_k=(d['price_yen']-500)/1000
            g,s,k,disc,bus=(d[z] for z in ['gate_code','safety_code','parking_service_code','bus_discount_code','bus_frequency_code'])
            score= (cf['reform_constant']*reform
                -(cf['price_nonuser_per_1000']+cf['additional_user_price_per_1000']*user)*p['_price_sensitivity']*price_k
                +cf['gate'][g]+cf['safety'][s]+cf['parking_service'][k]+cf['bus_discount'][disc]+cf['bus_frequency'][bus]
                + user*cf['user_parking_service_bonus'][k]
                +(1-user)*(cf['nonuser_bus_discount_bonus'][disc]+cf['nonuser_bus_frequency_bonus'][bus])
                +reform*cf['info_reform_bonus'][inf]
                +(cf['info_gate_bonus'] if inf=='I2' and g>0 else 0)
                +(cf['info_safety_bonus'] if inf=='I5' and s>0 else 0)
                +reform*p['_attitude']
            )
            # respondent-task-alternative shock shared between choice and support.
            shock=RNG.gauss(0, cf['shared_shock_sd'])
            score+=shock
            latent.append(score)
        # Draw a stochastic policy choice from the softmax probabilities.
        m=max(latent); w=[math.exp((u-m)/cf['choice_temperature']) for u in latent]
        chosen=sample_weighted([d['alternative'] for d in alts],w)
        row_task={'respondent_id':rid,'choice_set_id':f'{rid}_T{task:02d}',
                  'task_id':task,'block':p['block'],'survey_version':p['survey_version'],
                  'info_group':inf,'selected_alternative':chosen}
        for d,score in zip(alts,latent):
            alt=d['alternative']; reform=d['is_reform']; g,s,k,disc,bus=(d[z] for z in ['gate_code','safety_code','parking_service_code','bus_discount_code','bus_frequency_code'])
            if reform:
                # Ordered support: categories 1 oppose, 2 neutral, 3 support.
                support_latent=rf['intercept']+score*rf['shared_choice_score_weight']+RNG.gauss(0,rf['independent_shock_sd'])
                val=RNG.random()
                a=expit(rf['threshold_oppose']-support_latent);b=expit(rf['threshold_neutral']-support_latent)
                support=1 if val<a else (2 if val<b else 3)
                # Parking intention is NOT identical to policy support.
                index=(uf['base_user'] if user else uf['base_nonuser'])
                index+= -(uf['user_price_penalty'] if user else uf['nonuser_price_penalty'])*((d['price_yen']-500)/1000)
                index+=uf['safety_bonus'][s]+uf['parking_bonus'][k]-uf['bus_discount_penalty'][disc]-uf['bus_frequency_penalty'][bus]
                index+=p['_parking_propensity']+RNG.gauss(0,0.35)
                prob=expit(index)
                uncertain_prob=clip(uf['uncertainty_base']+uf['uncertainty_center_bonus']*(1-2*abs(prob-0.5)),0.01,0.40)
                intention=3 if RNG.random()<uncertain_prob else (1 if RNG.random()<prob else 2)
                row_task[f'{alt}_support']=support
                row_task[f'{alt}_parking_intention']=intention
            else:
                support='';intention=''
            out={'respondent_id':rid,'choice_set_id':f'{rid}_T{task:02d}',
                 'survey_version':p['survey_version'],'info_group':inf,'block':p['block'],
                 'task_id':task,'alternative':alt,'chosen':int(alt==chosen),
                 'support':support,'parking_intention':intention}
            choice_rows.append(out)
        task_rows.append(row_task)

# Public respondent file intentionally excludes latent generator-only variables.
public_people=[{k:v for k,v in p.items() if not k.startswith('_')} for p in people]
csv_write('respondents.csv',public_people)
csv_write('answers_per_task.csv',task_rows)
csv_write('responses_long.csv',choice_rows)

lookup={(d['task_id'],d['alternative']):d for d in design}
pmap={p['respondent_id']:p for p in public_people}
analysis=[]
for row in choice_rows:
    p=pmap[row['respondent_id']];d=lookup[(row['task_id'],row['alternative'])]
    dg,sa,ps,bd,bf=(d[z] for z in ['gate_code','safety_code','parking_service_code','bus_discount_code','bus_frequency_code'])
    reform=d['is_reform'];u=p['parking_user'];ig=p['info_group']
    a={**row,
       'price_yen':d['price_yen'],'price_1000':d['price_yen']/1000.0,
       'gate_code':dg,'safety_code':sa,'parking_service_code':ps,'bus_discount_code':bd,'bus_frequency_code':bf,'is_reform':reform,
       'gate_ic':int(dg==1),'gate_lpr':int(dg==2),
       'safety_basic':int(sa==1),'safety_plus':int(sa==2),
       'parking_live':int(ps==1),'parking_reserve':int(ps==2),
       'bus_disc50':int(bd==1),'bus_disc150':int(bd==2),
       'bus_plus1':int(bf==1),'bus_plus2':int(bf==2),
       'parking_user':u,'is_student':int(p['segment']=='学生'),'q1_role':p['q1_role'],
       'price_parking_user':(d['price_yen']/1000.0)*u,
       'parking_live_user':int(ps==1)*u,'parking_reserve_user':int(ps==2)*u,
       'bus_disc50_nonuser':int(bd==1)*(1-u),'bus_disc150_nonuser':int(bd==2)*(1-u),
       'bus_plus1_nonuser':int(bf==1)*(1-u),'bus_plus2_nonuser':int(bf==2)*(1-u),
       'reform_info2':reform*int(ig=='I2'),'reform_info3':reform*int(ig=='I3'),
       'reform_info4':reform*int(ig=='I4'),'reform_info5':reform*int(ig=='I5'),
       }
    analysis.append(a)
csv_write('analysis_long.csv',analysis)

# Automated content-level checks.
assert len(public_people)==1000 and len(task_rows)==6000 and len(choice_rows)==18000 and len(analysis)==18000
assert collections.Counter(p['segment'] for p in people)=={'学生':800,'教職員':200}
assert collections.Counter(p['parking_user'] for p in people)=={0:600,1:400}
assert all(collections.Counter(r['info_group'] for r in people)[f'I{i}']==200 for i in range(1,6))
assert all(66<=n<=67 for n in collections.Counter(p['survey_version'] for p in people).values())
assert all(v==1 for v in collections.Counter(r['choice_set_id'] for r in choice_rows if r['chosen']==1).values())
assert all(v==3 for v in collections.Counter(r['choice_set_id'] for r in choice_rows).values())
assert sum(1 for r in choice_rows if r['support']!='')==12000
assert all((r['support']=='' and r['parking_intention']=='') == (r['alternative']=='SQ') for r in choice_rows)
check={
 'seed':CONFIG['seed'],'respondents':len(public_people),'task_answers':len(task_rows),
 'long_rows':len(choice_rows),'support_answers':sum(1 for r in choice_rows if r['support']!=''),
 'users':sum(p['parking_user'] for p in people),'students':sum(p['segment']=='学生' for p in people),
 'choice_counts':dict(collections.Counter(r['selected_alternative'] for r in task_rows)),
 'support_counts':dict(collections.Counter(str(r['support']) for r in choice_rows if r['support']!='')),
 'intention_counts':dict(collections.Counter(str(r['parking_intention']) for r in choice_rows if r['parking_intention']!='')),
 'form_counts':dict(sorted(collections.Counter(p['survey_version'] for p in people).items())),
 'info_counts':dict(sorted(collections.Counter(p['info_group'] for p in people).items()))
}
(BASE/'validation_summary.json').write_text(json.dumps(check,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(check,ensure_ascii=False,indent=2))
