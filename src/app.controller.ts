import { BadRequestException, Controller, Get, HttpStatus, Post, Query } from '@nestjs/common';
import { AppService, ChartResult, Dataset } from './app.service';
import { Page } from './utils/decorators/page.decorator';
import { Redirect } from './utils/filters/redirect.filter';

@Controller()
export class AppController {
	constructor(private readonly service: AppService) {}

	@Page()
	@Get('/')
	public index(): PageProps {
		return {};
	}

	@Page()
	@Get('/compare')
	public comparison(
		@Query('datasets') datasets: string,
		@Query('gene') gene: string | undefined,
		@Query('groupBy') groupBy: string | undefined,
		@Query('splitBy') splitBy: string | undefined
	): CompareProps {
		const knownDatasets = this.getDatasets();

		if (!datasets.split(',').every((dataset) => knownDatasets.some((ds) => ds.name === dataset))) {
			throw new BadRequestException('Unknown dataset requested');
		}

		const genes = this.getGenes(datasets);

		if (gene === undefined || groupBy === undefined || splitBy === undefined) {
			const defaultGene = genes[0];

			if (defaultGene === undefined) {
				throw new BadRequestException('No common genes between datasets');
			}

			if (gene === undefined) gene = defaultGene;

			if (groupBy === undefined && splitBy === undefined) {
				groupBy = 'Genotype';
				splitBy = 'CellType';
			} else if (groupBy === undefined) {
				groupBy = splitBy === 'Genotype' ? 'CellType' : 'Genotype';
			} else if (splitBy === undefined) {
				splitBy = groupBy === 'Genotype' ? 'CellType' : 'Genotype';
			}

			throw new Redirect(`/compare?datasets=${datasets}&gene=${gene}&groupBy=${groupBy}&splitBy=${splitBy}`, HttpStatus.SEE_OTHER);
		}

		if (!genes.includes(gene)) {
			throw new BadRequestException('Selected gene is not in all datasets');
		}

		this.service.preload(datasets.split(','), gene, groupBy, splitBy);

		const sizes = Object.fromEntries(knownDatasets.map((ds) => [ds.name, ds.size]));

		return { order: datasets.split(',').sort((a, b) => sizes[a] - sizes[b]), genes, gene };
	}

	@Get('/datasets')
	public getDatasets(): Dataset[] {
		return this.service.getDatasets();
	}

	@Get('/genes')
	public getGenes(@Query('datasets') datasets: string): string[] {
		return this.service.getGenes(datasets.split(','));
	}

	@Get('/plot')
	public async getPlot(
		@Query('dataset') dataset: string,
		@Query('gene') gene: string,
		@Query('groupBy') groupBy: string,
		@Query('splitBy') splitBy: string
	): Promise<ChartResult> {
		return this.service.render(dataset, gene, groupBy, splitBy);
	}

	@Get('/plots')
	public async getPlots(
		@Query('datasets') datasets: string,
		@Query('gene') gene: string,
		@Query('groupBy') groupBy: string,
		@Query('splitBy') splitBy: string
	): Promise<Record<string, ChartResult>> {
		return this.service.generate(datasets.split(','), gene, groupBy, splitBy);
	}

	@Post('/preload')
	public async preload(@Query('dataset') datasets: string): Promise<void> {
		const genes = this.getGenes(datasets),
			defaultGene = genes[0];

		if (defaultGene === undefined) {
			throw new BadRequestException('No common genes between datasets');
		}

		this.service.preload(datasets.split(','), defaultGene, 'Genotype', 'CellType');
	}
}

