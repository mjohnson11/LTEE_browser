

const assays = ['Ara+2_R', 'Ara+2_2K', 'Ara+2_15K', 'Ara-1_R', 'Ara-1_2K', 'Ara-1_15K'];

const layout = {
  corr_div_left: 25,
  corr_div_top: 50,
  corr_div_width: 1100,
  corr_div_height: 500
}


const corr_graph_overlays = [];
let tooltip;
let gene_descrip_text;

let corr_data;
let ORF_list;
let count_ORF_map;
let count_data;
let row_clicked;

//colors
const hover_color = '#FF4444';
const click_color = '#4444FF';
const base_color = '#FCF0FC';
const background_color = '#000A00';

// axes etc.
const corr_domain = [-0.4, 0.1];

const corr_y = d3.scaleLinear().range([200, 60]).domain(corr_domain);
const log_freq_y = d3.scaleLinear().range([410, 280]).domain([-6, 0]);
const tn_x_scales = [];
const corr_x_scales = [];


// UTILITY AXES FUNCTION //

function standard_axes(svg_obj, xscale, yscale, xlabel, ylabel, title='') {
  const x_ax = d3.axisBottom(xscale).ticks(5);
  const y_ax = d3.axisLeft(yscale);
  if (yscale.domain().length == 2) y_ax.ticks(5);
  const xmid = (xscale.range()[0]+xscale.range()[1])/2;

  const y_called = svg_obj.append('g')
    .attr('class', 'axis')
    .attr('transform', 'translate(' + xscale.range()[0] + ', 0)')
    .call(y_ax);
  const x_called = svg_obj.append('g')
    .attr('class', 'axis')
    .attr('transform', 'translate(0,' + yscale.range()[0] + ')')
    .call(x_ax);

  if (yscale.domain().length > 2) update_axis(y_called);
  
  svg_obj.append('text')
    .attr('class', 'axis_label')
    .attr('text-anchor', 'middle')
    .attr('x', xmid)
    .attr('y', yscale.range()[0]+40)
    .text(xlabel)
    .attr('fill', base_color);

  svg_obj.append('text')
    .attr('class', 'axis_label')
    .attr('text-anchor', 'middle')
    .attr('x', xscale.range()[0])
    .attr('y', yscale.range()[1]-10)
    .text(ylabel)
    .attr('fill', base_color);

  svg_obj.append('text')
    .attr('class', 'axis_label')
    .attr('text-anchor', 'middle')
    .attr('x', xmid)
    .attr('y', yscale.range()[1]-20)
    .text(title)
    .attr('fill', base_color);
}

function make_search_bar() {
  search_div = d3.select('#main_div').append('div')
    .attr('id', 'search_div');
  gene_search = search_div.append('input')
    .attr('id', "gene_searchbar")
    .attr('type', 'search')
    .attr('spellcheck', 'false')
    .attr('autocomplete', 'off');

  //type="search" dir="ltr" spellcheck=false autocorrect="off" autocomplete="off" autocapitalize="off"
  const autoCompleteJS = new autoComplete({
      placeHolder: "Search for gene names...",
      selector: '#gene_searchbar',
      data: {src: ORF_list},
      events: {
          input: {
              selection: (event) => {
                  const selection = event.detail.selection.value;
                  autoCompleteJS.input.value = '';
                  d3.selectAll('.corr_group').classed('clicked_group', false);
                  let clicked = d3.selectAll('.corr_group').filter((d) => d['ORF']==selection);
                  clicked.raise();
                  clicked.selectAll('.corr_point').attr('r', 4);
                  clicked.classed('clicked_group', true);
                  document.getElementById('gene_searchbar').blur(); //removes focus so the cursor leaves
              }
          }
      }
  });

}



function make_tn_line(row_index) {
  let p = '';
  let datum;
  for (let count_data_row of count_ORF_map[row_index]) {
    datum = count_data[count_data_row];
    for (let i=0; i<6; i++) {
      let a = assays[i];
      let xs = tn_x_scales[i];
      p += 'M'+String(xs(1))+','+String(log_freq_y(parseFloat(datum[a+'_t'+String(1)+'_log10clipped'])));
      for (let t=2; t<6; t++) {
        //console.log(i, t, a+'_t'+String(t)+'_log10clipped', datum[a+'_t'+String(t)+'_log10clipped']);
        p += 'L'+String(xs(t))+','+String(log_freq_y(parseFloat(datum[a+'_t'+String(t)+'_log10clipped'])));
      }
    }
  }
  return p;
}

function one_corr_graph(svg_obj, xcol, ycol, xscale, yscale) {
  svg_obj.selectAll('.corr_group').append('circle')
    .attr('class', 'corr_point')
    .attr('r', 2)
    .attr('fill', base_color)
    .attr('cx', (d) => xscale(parseFloat(d[xcol])))
    .attr('cy', (d) => yscale(parseFloat(d[ycol])))
    .attr('opacity', (d) => ((d[xcol]=='') || (d[ycol]=='')) ? 0 : 1);
}

function make_corr_graphs() {
  // adding an index column
  for (let i=0; i< corr_data.length; i++) {
    corr_data[i]['WOW_index'] = i;
  }
  //corr_data = corr_data.slice(0, 10);

  const corr_div = d3.select('#data_div').append('div')
    .attr('class', 'graph_holder')
    .attr('id', 'corr_graphs_div')
    .style('left', layout.corr_div_left)
    .style('top', layout.corr_div_top)
    .style('width', layout.corr_div_width)
    .style('height', layout.corr_div_height);


  const corr_svg_obj = corr_div.append('svg')
    .attr('class', 'WOW_svg')
    .style('left', 0)
    .style('top', 0)
    .attr('width', layout.corr_div_width)
    .attr('height', layout.corr_div_height);


  corr_svg_obj.selectAll('.corr_group')
    .data(corr_data)
    .enter()
    .append('g')
      .attr('class', 'corr_group')
      .on('mouseover', function(e, d) {
        d3.selectAll('.clicked_group').raise();
        d3.select(this).raise().classed('hovered_group', true);
        d3.select(this).selectAll('.corr_point').attr('r', 4);
        show_tooltip(e.x, e.y, '<h2>'+d['ORF']+'</h2>');
      })
      .on('mouseout', function(e, d) {
        d3.select(this).selectAll('.corr_point').attr('r', 2);
        d3.selectAll('.clicked_group').selectAll('.corr_point').attr('r', 4);
        d3.select(this).classed('hovered_group', false);
        hide_tooltip();
      })
      .on('click', function(e, d) {
        d3.selectAll('.corr_group').classed('clicked_group', false);
        d3.select(this).selectAll('.corr_point').attr('r', 4);
        d3.select(this).classed('clicked_group', true);
      });
    
  
  for (let i=0; i<6; i++) {
    console.log(i);
    let range = [30+i*170, (i+1)*170-20];
    tn_x_scales.push(d3.scaleLinear().range(range).domain([0,5]));
    corr_x_scales.push(d3.scaleLinear().range(range).domain(corr_domain));
    standard_axes(corr_svg_obj, tn_x_scales[i], log_freq_y, 'Gens', i==0 ? 'log freq': '', '');
    standard_axes(corr_svg_obj, corr_x_scales[i], corr_y, 'S seg 1', i==0 ? 'S seg 2': '', assays[i]);
    one_corr_graph(corr_svg_obj, assays[i]+'_seg1', assays[i]+'_seg2', corr_x_scales[i], corr_y);
  }
  corr_svg_obj.selectAll('.corr_group').append('path')
    .attr('class', (d, i) => (i % 500 == 0) ? 'lin_path_example lin_path' : 'lin_path')
    .attr('d', (d, i) => make_tn_line(i));

}

function setup_tooltip() {
  tooltip = d3.select('body').append('div')
    .attr('class', 'WOW_tooltip')
    .html('<h2>yeah</h2><p>uhhuh</p>');
}

function show_tooltip(x, y, text) {
  tooltip.style('left', x-158).style('display', 'block').html(text);
  tooltip.style('top', y-8-tooltip.node().offsetHeight);
}

function hide_tooltip() {
  tooltip.style('display', 'none');
}

function load_data() {
  console.log('loading tnfit data');
  d3.csv('../Data_use/Couce_fitness_data.csv')
    .then(function(raw_data) {
      corr_data = raw_data; // using a global dataframe here
      ORF_list = corr_data.map((d) => d['ORF']);
      // load count data
      d3.csv('../Data_use/Couce_freq_data.csv')
      .then(function(raw_data) {
        count_data = raw_data;
        count_ORF_map = corr_data.map((d) => []);
        for (let i=0; i<count_data.length; i++) {
          let ii = ORF_list.indexOf(count_data[i]['ORF']);
          if (ii > -1) {
            count_ORF_map[ii].push(i);
          }
        }
        make_corr_graphs();
        setup_tooltip();
        make_search_bar();
        d3.select('#loading_text').remove();
        d3.select('#data_div').style('display', 'block');
      }).catch(function(error) {
        console.log(error);
      });
    })
    .catch(function(error) {
      console.log(error);
    });
}