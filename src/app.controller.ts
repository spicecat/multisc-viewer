import { BadRequestException, Controller, Get, HttpStatus, Query } from '@nestjs/common';
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
		@Query('groupBy') groupBy: string,
		@Query('splitBy') splitBy: string
	): CompareProps {
		const genes = this.getGenes(datasets);

		if (gene === undefined) {
			const defaultGene = genes[0];

			if (defaultGene === undefined) {
				throw new BadRequestException('No common genes between datasets');
			}

			throw new Redirect(`/compare?datasets=${datasets}&gene=${defaultGene}&groupBy=${groupBy}&splitBy=${splitBy}`, HttpStatus.SEE_OTHER);
		}

		if (!genes.includes(gene)) {
			throw new BadRequestException('Selected gene is not in all datasets');
		}

		this.service.preload(datasets.split(','), gene, groupBy, splitBy);

		return { genes, gene };
	}

	@Get('/datasets')
	public getDatasets(): Dataset[] {
		return this.service.getDatasets();
	}

	@Get('/genes')
	public getGenes(@Query('datasets') datasets: string): string[] {
		return this.service.getGenes(datasets.split(','));
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
}

